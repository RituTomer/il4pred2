#!/usr/bin/env python3
"""
IL4pred2 - Predicts IL-4 inducing peptides for human and mouse hosts
Developed by Prof G. P. S. Raghava's group
Please cite: https://webs.iiitd.edu.in/raghava/il4pred2/
"""

import argparse
import os
import re
import sys
from pathlib import Path
import shutil
from datetime import datetime
import numpy as np
import pandas as pd
import joblib
from subprocess import call
import warnings
warnings.filterwarnings('ignore')

def print_banner():
    """Display the tool banner."""
    banner = """
##############################################################################
# IL4pred2: Predicting IL-4 inducing peptides for human and mouse hosts      #
# Developed by Prof G. P. S. Raghava's group                                #
# Please cite: https://webs.iiitd.edu.in/raghava/il4pred2/                  #
##############################################################################
"""
    print(banner)

def initialize_environment():
    """Set up base directory and temporary folder."""
    base_dir = Path(__file__).parent.resolve()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    tmp_dir = base_dir / "tmp" / f"il4pred_{timestamp}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return base_dir, tmp_dir

def validate_sequence_length(sequence, min_len=9, max_len=22):
    """Validate that sequence length is between min_len and max_len."""
    length = len(sequence)
    if length < min_len or length > max_len:
        raise ValueError(
            f"Invalid sequence length {length}. "
            f"Sequence length must be between {min_len} and {max_len} amino acids."
        )
    return sequence

def read_sequences(input_file):
    """Read sequences from input file (FASTA or plain text)."""
    with open(input_file) as f:
        content = f.read().strip()
    
    if not content:
        raise ValueError("Input file is empty")

    # Handle FASTA format
    if content.startswith('>'):
        records = []
        current_id = None
        current_seq = []
        
        for line in content.split('\n'):
            if line.startswith('>'):
                if current_id is not None:
                    seq = ''.join(current_seq)
                    records.append((current_id, validate_sequence_length(seq)))
                current_id = line.split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        
        if current_id is not None:
            seq = ''.join(current_seq)
            records.append((current_id, validate_sequence_length(seq)))
    else:
        # Handle plain text format
        records = []
        for i, seq in enumerate(content.split('\n')):
            if seq.strip():
                cleaned_seq = seq.strip()
                records.append(
                    (f">Seq_{i+1}", validate_sequence_length(cleaned_seq))
                )
    
    # Validate and clean sequences
    cleaned_records = []
    for seq_id, seq in records:
        cleaned_seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', seq.upper())
        if not cleaned_seq:
            raise ValueError(f"Invalid sequence found in {seq_id}")
        cleaned_records.append((seq_id, cleaned_seq))
    
    return pd.DataFrame(cleaned_records, columns=['Sequence_ID', 'Sequence'])

def validate_file(path):
    """Validate that a file exists and is readable."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    if not os.access(path, os.R_OK):
        raise PermissionError(f"Cannot read file: {path}")
    return path

def get_model_path(host, job, base_dir):
    """Get the correct model path based on host and job type."""
    model_mapping = {
        (1, 1): "human_mlp_main.pkl",
        (1, 2): "human_et_alternate1_dpc.pkl",
        (1, 3): "human_et_alternate_2.pkl",
        (2, 1): "mouse_mlp_main.pkl",
        (2, 2): "mouse_lr_alternate_1.pkl",
        (2, 3): "mouse_rf_alternate_2.pkl"
    }
    try:
        model_file = base_dir / "model" / model_mapping[(host, job)]
        return validate_file(model_file)
    except KeyError:
        raise ValueError(f"Invalid host/job combination: host={host}, job={job}")

def generate_features(input_file, output_file, job, host, base_dir):
    """Generate features using external script."""
    input_path = validate_file(input_file)
    output_path = Path(output_file)
    feature_script = validate_file(base_dir / "pfeature_standalone" / "pfeature_comp.py")
    
    # Determine feature type based on host and job
    feature_types = {
        (1, 1): ("ALLCOMP", "human_main.csv"),
        (1, 2): ("DPC", None),
        (1, 3): ("ALLCOMP", "human_alt2.csv"),
        (2, 1): ("ALLCOMP", "mouse_main.csv"),
        (2, 2): ("ALLCOMP", "mouse_alt1.csv"),
        (2, 3): ("AAC", None)
    }
    
    feature_type, feature_file = feature_types.get((host, job), ("AAC", None))
    
    # Build feature command
    feature_cmd = [
        "python3", str(feature_script),
        "-i", str(input_path),
        "-o", str(output_path),
        "-j", feature_type
    ]
    
    # Execute feature generation
    original_dir = os.getcwd()
    try:
        os.chdir(feature_script.parent)
        if call(feature_cmd) != 0:
            raise RuntimeError("Feature generation failed")
        
        # Load and filter features if needed
        features = pd.read_csv(output_path)
        
        if feature_file:
            selected_features = pd.read_csv(
                base_dir / "Feature" / feature_file, 
                header=None
            )[0].tolist()
            try:
                features = features[selected_features]
            except KeyError as e:
                missing = set(selected_features) - set(features.columns)
                raise ValueError(f"Missing required features: {missing}")
        
        features.to_csv(output_path, index=False)
        return features.columns.tolist()
        
    finally:
        os.chdir(original_dir)

def run_prediction(features_file, model_file):
    """Run model prediction on generated features."""
    model = joblib.load(validate_file(model_file))
    features = pd.read_csv(features_file)
    
    # Ensure correct shape for prediction
    if len(features.shape) == 1:
        features = features.values.reshape(1, -1)
    else:
        features = features.values
    
    return model.predict_proba(features)[:, 1]

def save_results(sequences, scores, threshold, output_file, display_all=False):
    """Save results with predictions."""
    results = sequences.copy()
    results['Score'] = scores
    results['Prediction'] = np.where(
        scores >= threshold, 
        "IL4-inducing", 
        "IL4 non-inducing"
    )
    
    if not display_all:
        results = results[results['Prediction'] == "IL4-inducing"]
    
    results.to_csv(output_file, index=False)
    return results

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Predict IL-4 inducing peptides (9-22 amino acids)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input", 
        required=True,
        help="Input file with sequences (FASTA or one per line)"
    )
    parser.add_argument(
        "-o", "--output", 
        default="il4pred_results.csv",
        help="Output file for results"
    )
    parser.add_argument(
        "-s", "--host", 
        type=int, 
        choices=[1, 2], 
        default=1,
        help="Host: 1=Human, 2=Mouse"
    )
    parser.add_argument(
        "-j", "--job", 
        type=int, 
        choices=[1, 2, 3], 
        default=1,
        help="Prediction method: 1=Main, 2=Alt1, 3=Alt2"
    )
    parser.add_argument(
        "-t", "--threshold", 
        type=float, 
        default=0.5,
        help="Classification threshold (0-1)"
    )
    parser.add_argument(
        "-d", "--display", 
        type=int, 
        choices=[1, 2], 
        default=1,
        help="Display: 1=IL4 inducers only, 2=All peptides"
    )
    
    args = parser.parse_args()
    print_banner()
    
    # Initialize environment
    base_dir, tmp_dir = None, None
    try:
        base_dir, tmp_dir = initialize_environment()
        
        # Step 1: Read and validate input sequences
        print("\n[1/4] Reading and validating input sequences...")
        sequences = read_sequences(args.input)
        
        # Save validated sequences for feature generation
        seq_file = tmp_dir / "validated_sequences.fasta"
        sequences.to_csv(seq_file, sep='\n', header=False, index=False)
        
        # Step 2: Generate features
        print("[2/4] Generating features...")
        feature_file = tmp_dir / "features.csv"
        generate_features(seq_file, feature_file, args.job, args.host, base_dir)
        
        # Step 3: Run prediction
        print("[3/4] Running prediction...")
        model_file = get_model_path(args.host, args.job, base_dir)
        scores = run_prediction(feature_file, model_file)
        
        # Step 4: Save results
        print("[4/4] Saving results...")
        results = save_results(
            sequences, 
            scores, 
            args.threshold, 
            args.output, 
            display_all=(args.display == 2)
        )
        
        print(f"\nPrediction completed successfully!")
        print(f"Results saved to: {args.output}")
        print(f"IL4-inducing peptides found: {sum(results['Prediction'] == 'IL4-inducing')}")
        
    except Exception as e:
        print(f"\nERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)
    finally:
        # Always clean up temporary directory
        if tmp_dir and tmp_dir.exists():
            try:
                shutil.rmtree(tmp_dir, ignore_errors=True)
            except Exception as e:
                print(f"Warning: Failed to remove temporary directory {tmp_dir}: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()