# 🧬 תוכנית אימון מודל קלאסיפיקציה לסוגי גנים

## 📋 סיכום המשימה

**מטרה:** סיווג רצפי DNA ל-10 קטגוריות של סוגי גנים

| פרמטר | ערך |
|-------|-----|
| Input | NucleotideSequence (רצף DNA) |
| Output | GeneType (10 קטגוריות) |
| Train size | ~22,000 דוגמאות |
| אורך רצף | 2-1000 bp |

---

## 🔧 שלב 1: הכנת סביבת העבודה

### התקנת ספריות נדרשות

```bash
pip install pandas numpy scikit-learn tensorflow torch
pip install imbalanced-learn  # for SMOTE
pip install matplotlib seaborn  # for visualization
```

### מבנה תיקיות מומלץ

```
project/
├── data/
│   ├── train.csv
│   ├── validation.csv
│   └── test.csv
├── models/
│   └── saved_models/
├── notebooks/
│   └── exploration.ipynb
├── src/
│   ├── data_preprocessing.py
│   ├── feature_engineering.py
│   ├── model.py
│   └── train.py
└── results/
    └── metrics/
```

---

## 🧹 שלב 2: עיבוד מוקדם (Preprocessing)

### 2.1 טעינת נתונים

```python
import pandas as pd

# טעינה
train_df = pd.read_csv('train.csv')
val_df = pd.read_csv('validation.csv')
test_df = pd.read_csv('test.csv')

# שמירת רק עמודות רלוונטיות
columns_to_keep = ['NucleotideSequence', 'GeneType']
train_df = train_df[columns_to_keep]
```

### 2.2 ניקוי רצפים

```python
def clean_sequence(seq):
    """הסרת סימני < > וניקוי"""
    seq = seq.replace('<', '').replace('>', '')
    seq = seq.upper()
    # שמירה רק על ATGC
    seq = ''.join([c for c in seq if c in 'ATGC'])
    return seq

train_df['sequence'] = train_df['NucleotideSequence'].apply(clean_sequence)
```

### 2.3 סינון רצפים בעייתיים

```python
# הסרת רצפים קצרים מדי
MIN_LENGTH = 20
train_df = train_df[train_df['sequence'].str.len() >= MIN_LENGTH]

# הסרת כפילויות
train_df = train_df.drop_duplicates(subset=['sequence'])
```

### 2.4 טיפול בקטגוריות נדירות

```python
# אופציה 1: איחוד קטגוריות נדירות
rare_classes = ['scRNA']  # רק 3 דוגמאות!
train_df['GeneType'] = train_df['GeneType'].replace(rare_classes, 'OTHER')

# אופציה 2: הסרת קטגוריות נדירות
# train_df = train_df[~train_df['GeneType'].isin(rare_classes)]
```

---

## 🔢 שלב 3: הנדסת תכונות (Feature Engineering)

### 3.1 תכונות בסיסיות

```python
def extract_basic_features(seq):
    """חילוץ תכונות בסיסיות מרצף"""
    length = len(seq)
    
    # תדירות בסיסים
    a_freq = seq.count('A') / length if length > 0 else 0
    t_freq = seq.count('T') / length if length > 0 else 0
    g_freq = seq.count('G') / length if length > 0 else 0
    c_freq = seq.count('C') / length if length > 0 else 0
    
    # GC content
    gc_content = (g_freq + c_freq)
    
    return {
        'length': length,
        'a_freq': a_freq,
        't_freq': t_freq,
        'g_freq': g_freq,
        'c_freq': c_freq,
        'gc_content': gc_content
    }
```

### 3.2 K-mer Frequencies

```python
from itertools import product

def get_kmer_frequencies(seq, k=3):
    """חישוב תדירות k-mers"""
    kmers = [''.join(p) for p in product('ATGC', repeat=k)]
    kmer_counts = {kmer: 0 for kmer in kmers}
    
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in kmer_counts:
            kmer_counts[kmer] += 1
    
    # נרמול
    total = sum(kmer_counts.values())
    if total > 0:
        kmer_counts = {k: v/total for k, v in kmer_counts.items()}
    
    return kmer_counts
```

### 3.3 One-Hot Encoding לרצפים

```python
import numpy as np

def one_hot_encode(seq, max_length=500):
    """One-hot encoding לרצף DNA"""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    
    # Padding/Truncating
    if len(seq) > max_length:
        seq = seq[:max_length]
    
    encoded = np.zeros((max_length, 4))
    for i, char in enumerate(seq):
        if char in mapping:
            encoded[i, mapping[char]] = 1
    
    return encoded
```

---

## ⚖️ שלב 4: טיפול בחוסר איזון

### אופציה 1: Class Weights

```python
from sklearn.utils.class_weight import compute_class_weight

classes = train_df['GeneType'].unique()
weights = compute_class_weight('balanced', classes=classes, y=train_df['GeneType'])
class_weights = dict(zip(classes, weights))

# בזמן האימון
model.fit(X, y, class_weight=class_weights)
```

### אופציה 2: SMOTE (Oversampling)

```python
from imblearn.over_sampling import SMOTE

smote = SMOTE(random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_train, y_train)
```

### אופציה 3: Weighted Loss (Deep Learning)

```python
import torch.nn as nn

# חישוב משקלים
class_counts = train_df['GeneType'].value_counts()
weights = 1.0 / class_counts
weights = weights / weights.sum() * len(weights)

criterion = nn.CrossEntropyLoss(weight=torch.tensor(weights.values))
```

---

## 🤖 שלב 5: בחירת מודל

### אופציה A: Random Forest (Baseline)

```python
from sklearn.ensemble import RandomForestClassifier

# עם תכונות K-mer
model = RandomForestClassifier(
    n_estimators=200,
    max_depth=20,
    class_weight='balanced',
    random_state=42,
    n_jobs=-1
)

model.fit(X_train, y_train)
```

**יתרונות:** מהיר, interpretable, לא צריך GPU
**חסרונות:** לא לומד תבניות מורכבות ברצפים

### אופציה B: CNN 1D (מומלץ להתחלה)

```python
import tensorflow as tf
from tensorflow.keras import layers, models

def build_cnn_model(input_length=500, num_classes=10):
    model = models.Sequential([
        # Input: (batch, length, 4) - one-hot encoded
        layers.Input(shape=(input_length, 4)),
        
        # Conv blocks
        layers.Conv1D(64, 7, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(2),
        layers.Dropout(0.2),
        
        layers.Conv1D(128, 5, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.MaxPooling1D(2),
        layers.Dropout(0.2),
        
        layers.Conv1D(256, 3, activation='relu', padding='same'),
        layers.BatchNormalization(),
        layers.GlobalAveragePooling1D(),
        
        # Dense layers
        layers.Dense(128, activation='relu'),
        layers.Dropout(0.5),
        layers.Dense(num_classes, activation='softmax')
    ])
    
    return model

model = build_cnn_model()
model.compile(
    optimizer='adam',
    loss='sparse_categorical_crossentropy',
    metrics=['accuracy']
)
```

**יתרונות:** לומד תבניות מקומיות ברצף, יעיל
**חסרונות:** צריך GPU לאימון מהיר

### אופציה C: BiLSTM

```python
def build_lstm_model(input_length=500, num_classes=10):
    model = models.Sequential([
        layers.Input(shape=(input_length, 4)),
        
        layers.Bidirectional(layers.LSTM(64, return_sequences=True)),
        layers.Dropout(0.3),
        
        layers.Bidirectional(layers.LSTM(32)),
        layers.Dropout(0.3),
        
        layers.Dense(64, activation='relu'),
        layers.Dense(num_classes, activation='softmax')
    ])
    
    return model
```

**יתרונות:** לומד תלויות ארוכות
**חסרונות:** איטי יותר מ-CNN

### אופציה D: Transformer (מתקדם)

```python
def build_transformer_model(input_length=500, num_classes=10):
    inputs = layers.Input(shape=(input_length, 4))
    
    # Positional encoding
    x = layers.Dense(64)(inputs)
    
    # Transformer block
    attention = layers.MultiHeadAttention(num_heads=4, key_dim=64)(x, x)
    x = layers.Add()([x, attention])
    x = layers.LayerNormalization()(x)
    
    ff = layers.Dense(128, activation='relu')(x)
    ff = layers.Dense(64)(ff)
    x = layers.Add()([x, ff])
    x = layers.LayerNormalization()(x)
    
    # Output
    x = layers.GlobalAveragePooling1D()(x)
    x = layers.Dense(64, activation='relu')(x)
    outputs = layers.Dense(num_classes, activation='softmax')(x)
    
    return models.Model(inputs, outputs)
```

**יתרונות:** State-of-the-art לרצפים
**חסרונות:** צריך הרבה נתונים ומשאבים

---

## 📊 שלב 6: אימון והערכה

### 6.1 הגדרת Callbacks

```python
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

callbacks = [
    EarlyStopping(
        monitor='val_loss',
        patience=10,
        restore_best_weights=True
    ),
    ModelCheckpoint(
        'best_model.h5',
        monitor='val_loss',
        save_best_only=True
    ),
    ReduceLROnPlateau(
        monitor='val_loss',
        factor=0.5,
        patience=5
    )
]
```

### 6.2 אימון

```python
history = model.fit(
    X_train, y_train,
    validation_data=(X_val, y_val),
    epochs=100,
    batch_size=32,
    class_weight=class_weights,
    callbacks=callbacks
)
```

### 6.3 הערכה

```python
from sklearn.metrics import classification_report, confusion_matrix, f1_score

# ניבוי
y_pred = model.predict(X_test)
y_pred_classes = np.argmax(y_pred, axis=1)

# מדדים
print("F1 Score (Macro):", f1_score(y_test, y_pred_classes, average='macro'))
print("\nClassification Report:")
print(classification_report(y_test, y_pred_classes, target_names=class_names))

# Confusion Matrix
cm = confusion_matrix(y_test, y_pred_classes)
```

---

## 🎯 שלב 7: אופטימיזציה

### 7.1 Hyperparameter Tuning

```python
from sklearn.model_selection import GridSearchCV

# For Random Forest
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [10, 20, 30, None],
    'min_samples_split': [2, 5, 10]
}

grid_search = GridSearchCV(
    RandomForestClassifier(class_weight='balanced'),
    param_grid,
    cv=5,
    scoring='f1_macro',
    n_jobs=-1
)
```

### 7.2 Ensemble Methods

```python
from sklearn.ensemble import VotingClassifier

ensemble = VotingClassifier([
    ('rf', RandomForestClassifier(n_estimators=200, class_weight='balanced')),
    ('xgb', XGBClassifier(scale_pos_weight=...)),
    ('lgb', LGBMClassifier(class_weight='balanced'))
], voting='soft')
```

---

## 📅 לוח זמנים מומלץ

| שלב | משימה | זמן משוער |
|-----|-------|-----------|
| 1 | הכנת סביבה | 30 דקות |
| 2 | עיבוד מוקדם | 1-2 שעות |
| 3 | הנדסת תכונות | 2-3 שעות |
| 4 | Baseline (Random Forest) | 1 שעה |
| 5 | CNN Model | 2-4 שעות |
| 6 | הערכה ושיפור | 2-3 שעות |
| 7 | Fine-tuning | 2-4 שעות |

**סה"כ:** ~1-2 ימי עבודה

---

## ⚠️ נקודות חשובות לזכור

### ❌ מה לא לעשות:
1. **לא להשתמש ב-Symbol/Description** - דליפת מידע!
2. **לא להשתמש ב-Accuracy** - מטעה בגלל חוסר איזון
3. **לא להתעלם מהחפיפה בין Splits** - יוצר תוצאות מוטות
4. **לא להתעלם מקטגוריות נדירות** - scRNA עם 3 דוגמאות

### ✅ מה כן לעשות:
1. **להשתמש ב-F1 Macro** כמדד עיקרי
2. **לטפל בחוסר איזון** - class weights / SMOTE
3. **ליצור splits חדשים** ללא חפיפה
4. **לאחד קטגוריות נדירות** (scRNA → OTHER)
5. **לשמור מודלים** עם checkpoints
6. **לתעד ניסויים** עם MLflow/WandB

---

## 📈 ציפיות לביצועים

| גישה | F1-Macro צפוי |
|------|---------------|
| Baseline (Majority) | ~0.10 |
| Random Forest + K-mers | 0.40-0.55 |
| CNN 1D | 0.60-0.75 |
| BiLSTM | 0.65-0.80 |
| Transformer | 0.70-0.85 |
| Ensemble | 0.75-0.90 |

---

## 🚀 המלצה: התחל כאן!

```python
# Quick Start - התחלה מהירה
# 1. Random Forest עם K-mers כ-baseline
# 2. CNN פשוט אם יש GPU
# 3. שפר בהדרגה
```

---

*תוכנית זו נוצרה ב: 2024-12-11*

