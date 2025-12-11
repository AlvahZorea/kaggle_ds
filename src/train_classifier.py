#!/usr/bin/env python3
"""
Gene Type Classification Model Training Script
סקריפט אימון מודל קלאסיפיקציה לסוגי גנים

Usage:
    python train_classifier.py --model [rf|cnn|lstm] --epochs 50

Requirements:
    pip install pandas numpy scikit-learn tensorflow imbalanced-learn
"""

import argparse
import os
import numpy as np
import pandas as pd
from collections import Counter
from itertools import product

# Sklearn
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, f1_score
from sklearn.utils.class_weight import compute_class_weight

# For handling imbalance
from imblearn.over_sampling import SMOTE

import warnings
warnings.filterwarnings('ignore')


# ============================================================
# Data Loading & Preprocessing
# ============================================================

def load_data(data_dir='./'):
    """טעינת נתונים"""
    print("📂 Loading data...")
    
    train_df = pd.read_csv(os.path.join(data_dir, 'train.csv'))
    val_df = pd.read_csv(os.path.join(data_dir, 'validation.csv'))
    test_df = pd.read_csv(os.path.join(data_dir, 'test.csv'))
    
    print(f"   Train: {len(train_df)} samples")
    print(f"   Val: {len(val_df)} samples")
    print(f"   Test: {len(test_df)} samples")
    
    return train_df, val_df, test_df


def clean_sequence(seq):
    """ניקוי רצף DNA"""
    if pd.isna(seq):
        return ''
    seq = str(seq).replace('<', '').replace('>', '').upper()
    seq = ''.join([c for c in seq if c in 'ATGC'])
    return seq


def preprocess_data(df, min_length=20):
    """עיבוד מוקדם"""
    df = df.copy()
    
    # ניקוי רצפים
    df['sequence'] = df['NucleotideSequence'].apply(clean_sequence)
    
    # סינון רצפים קצרים
    df = df[df['sequence'].str.len() >= min_length]
    
    # איחוד קטגוריות נדירות
    rare_classes = ['scRNA']
    df['GeneType'] = df['GeneType'].replace(rare_classes, 'OTHER')
    
    return df


def remove_duplicates_between_splits(train_df, val_df, test_df):
    """הסרת חפיפות בין הסטים"""
    print("🧹 Removing duplicates between splits...")
    
    train_seqs = set(train_df['sequence'])
    val_seqs = set(val_df['sequence'])
    test_seqs = set(test_df['sequence'])
    
    # הסרת רצפים מ-val/test שמופיעים ב-train
    val_df = val_df[~val_df['sequence'].isin(train_seqs)]
    test_df = test_df[~test_df['sequence'].isin(train_seqs)]
    
    # הסרת רצפים מ-test שמופיעים ב-val
    test_df = test_df[~test_df['sequence'].isin(val_seqs)]
    
    print(f"   After cleaning - Train: {len(train_df)}, Val: {len(val_df)}, Test: {len(test_df)}")
    
    return train_df, val_df, test_df


# ============================================================
# Feature Engineering
# ============================================================

def extract_basic_features(seq):
    """חילוץ תכונות בסיסיות"""
    length = len(seq)
    if length == 0:
        return [0] * 6
    
    a_freq = seq.count('A') / length
    t_freq = seq.count('T') / length
    g_freq = seq.count('G') / length
    c_freq = seq.count('C') / length
    gc_content = g_freq + c_freq
    
    return [length, a_freq, t_freq, g_freq, c_freq, gc_content]


def get_kmer_frequencies(seq, k=3):
    """חישוב תדירות k-mers"""
    kmers = [''.join(p) for p in product('ATGC', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in kmers}
    
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
    
    total = sum(kmer_counts.values())
    if total > 0:
        return [kmer_counts[kmer]/total for kmer in kmers]
    return [0] * len(kmers)


def extract_features(sequences, use_kmers=True, kmer_k=3):
    """חילוץ כל התכונות"""
    print("🔧 Extracting features...")
    
    features = []
    for seq in sequences:
        feat = extract_basic_features(seq)
        if use_kmers:
            feat.extend(get_kmer_frequencies(seq, k=kmer_k))
        features.append(feat)
    
    return np.array(features)


def one_hot_encode_sequences(sequences, max_length=500):
    """One-hot encoding לרצפים (עבור CNN/LSTM)"""
    print("🔢 One-hot encoding sequences...")
    
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    encoded = np.zeros((len(sequences), max_length, 4))
    
    for i, seq in enumerate(sequences):
        for j, char in enumerate(seq[:max_length]):
            if char in mapping:
                encoded[i, j, mapping[char]] = 1
    
    return encoded


# ============================================================
# Models
# ============================================================

def train_random_forest(X_train, y_train, X_val, y_val, class_weights=None):
    """אימון Random Forest"""
    print("🌲 Training Random Forest...")
    
    model = RandomForestClassifier(
        n_estimators=200,
        max_depth=20,
        min_samples_split=5,
        class_weight='balanced',
        random_state=42,
        n_jobs=-1
    )
    
    model.fit(X_train, y_train)
    
    # הערכה על validation
    y_pred = model.predict(X_val)
    f1 = f1_score(y_val, y_pred, average='macro')
    print(f"   Validation F1-Macro: {f1:.4f}")
    
    return model


def build_cnn_model(input_length=500, num_classes=9):
    """בניית מודל CNN"""
    try:
        import tensorflow as tf
        from tensorflow.keras import layers, models
    except ImportError:
        print("❌ TensorFlow not installed. Run: pip install tensorflow")
        return None
    
    model = models.Sequential([
        layers.Input(shape=(input_length, 4)),
        
        # Conv Block 1
        layers.Conv1D(64, 7, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(2),
        layers.Dropout(0.2),
        
        # Conv Block 2
        layers.Conv1D(128, 5, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(2),
        layers.Dropout(0.2),
        
        # Conv Block 3
        layers.Conv1D(256, 3, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.GlobalAveragePooling1D(),
        
        # Dense
        layers.Dense(128, activation='relu'),
        layers.Dropout(0.5),
        layers.Dense(num_classes, activation='softmax')
    ])
    
    model.compile(
        optimizer='adam',
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    return model


def train_cnn(X_train, y_train, X_val, y_val, class_weights, epochs=50):
    """אימון מודל CNN"""
    try:
        from tensorflow.keras.callbacks import EarlyStopping, ReduceLROnPlateau
    except ImportError:
        print("❌ TensorFlow not installed.")
        return None
    
    print("🧠 Training CNN...")
    
    num_classes = len(np.unique(y_train))
    model = build_cnn_model(input_length=X_train.shape[1], num_classes=num_classes)
    
    if model is None:
        return None
    
    callbacks = [
        EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True),
        ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=5)
    ]
    
    history = model.fit(
        X_train, y_train,
        validation_data=(X_val, y_val),
        epochs=epochs,
        batch_size=32,
        class_weight=class_weights,
        callbacks=callbacks,
        verbose=1
    )
    
    return model


def build_lstm_model(input_length=500, num_classes=9):
    """בניית מודל LSTM"""
    try:
        from tensorflow.keras import layers, models
    except ImportError:
        print("❌ TensorFlow not installed.")
        return None
    
    model = models.Sequential([
        layers.Input(shape=(input_length, 4)),
        layers.Bidirectional(layers.LSTM(64, return_sequences=True)),
        layers.Dropout(0.3),
        layers.Bidirectional(layers.LSTM(32)),
        layers.Dropout(0.3),
        layers.Dense(64, activation='relu'),
        layers.Dense(num_classes, activation='softmax')
    ])
    
    model.compile(
        optimizer='adam',
        loss='sparse_categorical_crossentropy',
        metrics=['accuracy']
    )
    
    return model


# ============================================================
# Evaluation
# ============================================================

def evaluate_model(model, X_test, y_test, label_encoder, model_type='rf'):
    """הערכת מודל"""
    print("\n📊 Evaluating model...")
    
    if model_type == 'rf':
        y_pred = model.predict(X_test)
    else:
        y_pred_proba = model.predict(X_test)
        y_pred = np.argmax(y_pred_proba, axis=1)
    
    # מדדים
    f1_macro = f1_score(y_test, y_pred, average='macro')
    f1_weighted = f1_score(y_test, y_pred, average='weighted')
    
    print(f"\n{'='*50}")
    print(f"F1-Score (Macro): {f1_macro:.4f}")
    print(f"F1-Score (Weighted): {f1_weighted:.4f}")
    print(f"{'='*50}")
    
    # Classification report
    class_names = label_encoder.classes_
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred, target_names=class_names))
    
    # Confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    print("\nConfusion Matrix:")
    print(cm)
    
    return f1_macro


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='Train Gene Type Classifier')
    parser.add_argument('--model', type=str, default='rf', choices=['rf', 'cnn', 'lstm'],
                        help='Model type: rf (Random Forest), cnn, lstm')
    parser.add_argument('--epochs', type=int, default=50, help='Number of epochs for CNN/LSTM')
    parser.add_argument('--max-length', type=int, default=500, help='Max sequence length')
    parser.add_argument('--data-dir', type=str, default='./', help='Data directory')
    args = parser.parse_args()
    
    print("=" * 60)
    print("🧬 Gene Type Classification Training")
    print("=" * 60)
    
    # 1. Load data
    train_df, val_df, test_df = load_data(args.data_dir)
    
    # 2. Preprocess
    print("\n🧹 Preprocessing...")
    train_df = preprocess_data(train_df)
    val_df = preprocess_data(val_df)
    test_df = preprocess_data(test_df)
    
    # 3. Remove duplicates between splits
    train_df, val_df, test_df = remove_duplicates_between_splits(train_df, val_df, test_df)
    
    # 4. Encode labels
    label_encoder = LabelEncoder()
    y_train = label_encoder.fit_transform(train_df['GeneType'])
    y_val = label_encoder.transform(val_df['GeneType'])
    y_test = label_encoder.transform(test_df['GeneType'])
    
    print(f"\nClasses: {label_encoder.classes_}")
    print(f"Class distribution (train): {Counter(y_train)}")
    
    # 5. Compute class weights
    class_weights = compute_class_weight('balanced', classes=np.unique(y_train), y=y_train)
    class_weight_dict = dict(zip(np.unique(y_train), class_weights))
    print(f"Class weights: {class_weight_dict}")
    
    # 6. Extract features / encode sequences
    if args.model == 'rf':
        print("\n" + "="*50)
        print("Using K-mer features for Random Forest")
        print("="*50)
        
        X_train = extract_features(train_df['sequence'].values, use_kmers=True, kmer_k=3)
        X_val = extract_features(val_df['sequence'].values, use_kmers=True, kmer_k=3)
        X_test = extract_features(test_df['sequence'].values, use_kmers=True, kmer_k=3)
        
        model = train_random_forest(X_train, y_train, X_val, y_val)
        evaluate_model(model, X_test, y_test, label_encoder, model_type='rf')
        
    else:  # cnn or lstm
        print("\n" + "="*50)
        print(f"Using One-Hot Encoding for {args.model.upper()}")
        print("="*50)
        
        X_train = one_hot_encode_sequences(train_df['sequence'].values, max_length=args.max_length)
        X_val = one_hot_encode_sequences(val_df['sequence'].values, max_length=args.max_length)
        X_test = one_hot_encode_sequences(test_df['sequence'].values, max_length=args.max_length)
        
        if args.model == 'cnn':
            model = train_cnn(X_train, y_train, X_val, y_val, class_weight_dict, epochs=args.epochs)
        else:  # lstm
            model = build_lstm_model(input_length=args.max_length, num_classes=len(label_encoder.classes_))
            if model:
                from tensorflow.keras.callbacks import EarlyStopping
                model.fit(X_train, y_train, validation_data=(X_val, y_val),
                         epochs=args.epochs, batch_size=32, class_weight=class_weight_dict,
                         callbacks=[EarlyStopping(patience=10, restore_best_weights=True)])
        
        if model:
            evaluate_model(model, X_test, y_test, label_encoder, model_type=args.model)
    
    print("\n✅ Training complete!")


if __name__ == "__main__":
    main()

