#!/usr/bin/env python3
"""
סקריפט לניתוח נתונים וליצירת דוח מסכם של המשתנים והתיוגים
"""

import pandas as pd
import numpy as np
from collections import Counter
import os

def analyze_dataset(file_path, dataset_name):
    """מנתח קובץ נתונים ומחזיר סטטיסטיקות"""
    print(f"\n{'='*80}")
    print(f"ניתוח קובץ: {dataset_name}")
    print(f"{'='*80}")
    
    # קריאת הקובץ
    df = pd.read_csv(file_path)
    
    print(f"\n📊 מידע כללי:")
    print(f"  • מספר שורות: {len(df):,}")
    print(f"  • מספר עמודות: {len(df.columns)}")
    
    print(f"\n📋 רשימת המשתנים (עמודות):")
    for i, col in enumerate(df.columns, 1):
        print(f"  {i}. {col}")
    
    print(f"\n🔍 סוגי הנתונים:")
    for col in df.columns:
        dtype = df[col].dtype
        null_count = df[col].isnull().sum()
        null_pct = (null_count / len(df)) * 100
        print(f"  • {col}: {dtype} (ערכים חסרים: {null_count} ({null_pct:.2f}%))")
    
    # ניתוח התיוגים (GeneType)
    if 'GeneType' in df.columns:
        print(f"\n🏷️  ניתוח התיוגים (GeneType):")
        label_counts = df['GeneType'].value_counts()
        print(f"  • מספר קטגוריות שונות: {len(label_counts)}")
        print(f"\n  התפלגות התיוגים:")
        for label, count in label_counts.items():
            pct = (count / len(df)) * 100
            print(f"    - {label}: {count:,} ({pct:.2f}%)")
    
    # ניתוח GeneGroupMethod
    if 'GeneGroupMethod' in df.columns:
        print(f"\n🔬 ניתוח GeneGroupMethod:")
        method_counts = df['GeneGroupMethod'].value_counts()
        for method, count in method_counts.items():
            pct = (count / len(df)) * 100
            print(f"  • {method}: {count:,} ({pct:.2f}%)")
    
    # ניתוח אורך רצפי הנוקלאוטידים
    if 'NucleotideSequence' in df.columns:
        print(f"\n🧬 ניתוח רצפי נוקלאוטידים:")
        # הסרת הסימנים < > מהרצפים
        sequences = df['NucleotideSequence'].str.replace('<', '').str.replace('>', '')
        lengths = sequences.str.len()
        print(f"  • אורך ממוצע: {lengths.mean():.2f} נוקלאוטידים")
        print(f"  • אורך מינימלי: {lengths.min()} נוקלאוטידים")
        print(f"  • אורך מקסימלי: {lengths.max()} נוקלאוטידים")
        print(f"  • חציון: {lengths.median():.2f} נוקלאוטידים")
        print(f"  • סטיית תקן: {lengths.std():.2f} נוקלאוטידים")
    
    # ניתוח NCBIGeneID
    if 'NCBIGeneID' in df.columns:
        print(f"\n🆔 ניתוח NCBIGeneID:")
        unique_ids = df['NCBIGeneID'].nunique()
        print(f"  • מספר ID ייחודיים: {unique_ids:,}")
        print(f"  • מספר כפילויות: {len(df) - unique_ids:,}")
    
    # ניתוח Symbol
    if 'Symbol' in df.columns:
        print(f"\n📝 ניתוח Symbol:")
        unique_symbols = df['Symbol'].nunique()
        print(f"  • מספר סמלים ייחודיים: {unique_symbols:,}")
        print(f"  • מספר כפילויות: {len(df) - unique_symbols:,}")
    
    return df

def create_summary_report():
    """יוצר דוח מסכם של כל קבצי הנתונים"""
    print("\n" + "="*80)
    print("דוח מסכם - ניתוח נתונים")
    print("="*80)
    
    datasets = {
        'train.csv': 'אימון (Train)',
        'validation.csv': 'ולידציה (Validation)',
        'test.csv': 'בדיקה (Test)'
    }
    
    all_data = {}
    
    for file_name, display_name in datasets.items():
        if os.path.exists(file_name):
            df = analyze_dataset(file_name, display_name)
            all_data[display_name] = df
        else:
            print(f"\n⚠️  קובץ {file_name} לא נמצא")
    
    # השוואה בין הקבצים
    if len(all_data) > 1:
        print(f"\n{'='*80}")
        print("השוואה בין קבצי הנתונים")
        print(f"{'='*80}")
        
        comparison_data = []
        for name, df in all_data.items():
            comparison_data.append({
                'קובץ': name,
                'מספר שורות': len(df),
                'מספר עמודות': len(df.columns)
            })
        
        comparison_df = pd.DataFrame(comparison_data)
        print("\n" + comparison_df.to_string(index=False))
        
        # השוואת התפלגות התיוגים
        if 'GeneType' in list(all_data.values())[0].columns:
            print(f"\n{'='*80}")
            print("השוואת התפלגות התיוגים בין הקבצים")
            print(f"{'='*80}")
            
            for name, df in all_data.items():
                print(f"\n{name}:")
                label_counts = df['GeneType'].value_counts()
                for label, count in label_counts.head(10).items():
                    pct = (count / len(df)) * 100
                    print(f"  {label}: {count:,} ({pct:.2f}%)")

if __name__ == "__main__":
    create_summary_report()

