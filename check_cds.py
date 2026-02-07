#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CDS Codon Usage Checker v1.0
Validate FASTA format CDS sequences and generate HTML reports

Author: Liu Haoqiu
Date: 2025-09-08
"""

import os
import sys
import subprocess
import importlib

def check_and_install_packages():
    """
    Check and automatically install required Python packages
    """
    required_packages = {
        'matplotlib': 'matplotlib',
        'numpy': 'numpy', 
        'pandas': 'pandas'
    }
    
    missing_packages = []
    
    print("Checking required Python packages...")
    
    for package_name, pip_name in required_packages.items():
        try:
            importlib.import_module(package_name)
            print(f"✅ {package_name} installed")
        except ImportError:
            print(f"❌ {package_name} not installed")
            missing_packages.append(pip_name)
        except Exception as e:
            print(f"⚠️ Exception while checking {package_name}: {e}")
            missing_packages.append(pip_name)
    
    if missing_packages:
        print(f"\nFound {len(missing_packages)} missing packages. Installing...")
        for package in missing_packages:
            print(f"Installing {package}...")
            try:
                subprocess.check_call([sys.executable, '-m', 'pip', 'install', package], 
                                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                print(f"✅ {package} installed successfully")
            except subprocess.CalledProcessError as e:
                print(f"❌ Failed to install {package}: {e}")
                print("Please install the missing packages manually and rerun the program")
                sys.exit(1)
        print("\nAll packages installed.")
    else:
        print("✅ All required packages are installed")
    
    print("-" * 50)

# Check dependencies before importing other packages
check_and_install_packages()

# Now safely import all packages
import os
import re
import time
import html
from datetime import datetime
from collections import defaultdict, Counter
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import json


class CDSChecker:
    """CDS sequence validator"""
    
    def __init__(self, input_file="./cds_from_genomic.fna"):
        self.input_file = input_file
        self.version = "v1.22"
        self.timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        # Statistical data
        self.total_records = 0
        self.valid_records = 0
        self.invalid_records = []
        self.non_atcg_records = []
        self.valid_lengths = []
        
        # Violation statistics
        self.violations = {
            'length': 0,
            'start_codon': 0,
            'stop_codon': 0
        }
        
        # Valid start and stop codons
        self.valid_start_codons = {'ATG'}
        self.valid_stop_codons = {'TGA', 'TAA', 'TAG'}
        self.valid_nucleotides = set('ATCG')
        
        # Synonymous codon mapping (including STOP group)
        self.codon_to_aa = {
            # Ala (A)
            'GCT': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
            # Arg (R)
            'CGT': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
            # Asn (N)
            'AAT': 'Asn (N)', 'AAC': 'Asn (N)',
            # Asp (D)
            'GAT': 'Asp (D)', 'GAC': 'Asp (D)',
            # Cys (C)
            'TGT': 'Cys (C)', 'TGC': 'Cys (C)',
            # Gln (Q)
            'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
            # Glu (E)
            'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
            # Gly (G)
            'GGT': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)',
            # His (H)
            'CAT': 'His (H)', 'CAC': 'His (H)',
            # Ile (I)
            'ATT': 'Ile (I)', 'ATC': 'Ile (I)', 'ATA': 'Ile (I)',
            # Leu (L)
            'CTT': 'Leu (L)', 'CTC': 'Leu (L)', 'CTA': 'Leu (L)', 'CTG': 'Leu (L)', 'TTA': 'Leu (L)', 'TTG': 'Leu (L)',
            # Lys (K)
            'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
            # Met (M)
            'ATG': 'Met (M)',
            # Phe (F)
            'TTT': 'Phe (F)', 'TTC': 'Phe (F)',
            # Pro (P)
            'CCT': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
            # Ser (S)
            'TCT': 'Ser (S)', 'TCC': 'Ser (S)', 'TCA': 'Ser (S)', 'TCG': 'Ser (S)', 'AGT': 'Ser (S)', 'AGC': 'Ser (S)',
            # Thr (T)
            'ACT': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
            # Trp (W)
            'TGG': 'Trp (W)',
            # Tyr (Y)
            'TAT': 'Tyr (Y)', 'TAC': 'Tyr (Y)',
            # Val (V)
            'GTT': 'Val (V)', 'GTC': 'Val (V)', 'GTA': 'Val (V)', 'GTG': 'Val (V)',
            # STOP
            'TAA': 'STOP', 'TGA': 'STOP', 'TAG': 'STOP'
        }
        
        # Standard 64 codon set
        self.standard_codons = set(self.codon_to_aa.keys())
        
        # Codon usage statistics
        self.participating_genes = 0
        self.excluded_atcg_bad_genes = 0
        self.total_codon_count = 0
        self.codon_counts = defaultdict(lambda: defaultdict(int))
        self.aa_totals = defaultdict(int)
        self.codon_freq = defaultdict(lambda: defaultdict(float))
        self.skipped_triplets = 0
    
    def parse_fasta(self):
        """
        Parse FASTA file, reading entire file content at once
        """
        try:
            with open(self.input_file, 'r', encoding='utf-8') as f:
                content = f.read()
        except FileNotFoundError:
            print(f"Error: File {self.input_file} not found")
            return []
        except Exception as e:
            print(f"Error: Exception occurred while reading file - {e}")
            return []
        
        if not content.strip():
            print("Warning: Input file is empty")
            return []
        
        # Split records by '>'
        records = []
        parts = content.split('>')
        
        for part in parts[1:]:  # skip the first empty part
            lines = part.strip().split('\n')
            if not lines:
                continue
                
            header = lines[0].strip()
            if not header:
                continue
                
            # Merge all sequence lines
            sequence_lines = lines[1:]
            sequence = ''.join(line.strip() for line in sequence_lines)
            
            # Clean sequence: convert to uppercase and remove whitespace
            sequence = ''.join(sequence.upper().split())
            
            if sequence:  # Only add when sequence is not empty
                records.append((header, sequence))
        
        return records
    
    def detect_non_atcg_characters(self, sequence):
        """
        Detect non-ATCG characters in sequence
        """
        sequence_set = set(sequence)
        non_atcg = sequence_set - self.valid_nucleotides
        return sorted(list(non_atcg)) if non_atcg else []
    
    def validate_cds(self, header, sequence):
        """
        Validate compliance of a single CDS sequence
        Returns: (is_valid, violations, warnings, start_codon, stop_codon)
        """
        violations = []
        warnings = []
        start_codon = ""
        stop_codon = ""
        
        seq_length = len(sequence)
        
        # Detect non-ATCG characters
        non_atcg_chars = self.detect_non_atcg_characters(sequence)
        if non_atcg_chars:
            warnings.append("Contains non-ATCG symbols")
            self.non_atcg_records.append({
                'gene': header,
                'length': seq_length,
                'non_atcg_chars': ', '.join(non_atcg_chars)
            })
        
        # Rule 1: Sequence length must be a multiple of 3
        length_valid = (seq_length % 3 == 0)
        if not length_valid:
            violations.append("len%3!=0")
            self.violations['length'] += 1
        
        # Rule 2: Start codon must be ATG
        if seq_length >= 3:
            start_codon = sequence[:3]
            if start_codon not in self.valid_start_codons:
                violations.append("start≠ATG")
                self.violations['start_codon'] += 1
        else:
            violations.append("start≠ATG")
            self.violations['start_codon'] += 1
        
        # Rule 3: Stop codon check (only when length is multiple of 3)
        if length_valid and seq_length >= 3:
            stop_codon = sequence[-3:]
            if stop_codon not in self.valid_stop_codons:
                violations.append("stop∉{TGA,TAA,TAG}")
                self.violations['stop_codon'] += 1
        elif not length_valid:
            # When length is not multiple of 3, skip stop codon check but record as empty
            pass
        
        is_valid = len(violations) == 0
        
        return is_valid, violations, warnings, start_codon, stop_codon
    
    def is_pure_atcg(self, sequence):
        """
        Check if sequence contains only ATCG characters
        """
        sequence_set = set(sequence)
        return sequence_set.issubset(self.valid_nucleotides)
    
    def analyze_codon_usage(self, sequence):
        """
        Analyze codon usage of a single pure ATCG compliant CDS sequence
        Note: This method only processes pure ATCG sequences, thus no skipped triplets
        """
        seq_length = len(sequence)
        
        # Split sequence to end with step size 3 (including stop codon)
        for i in range(0, seq_length, 3):
            if i + 3 <= seq_length:
                triplet = sequence[i:i+3]
                
                # Since input is pre-filtered to pure ATCG, all 3bp triplets should be standard codons
                if triplet in self.standard_codons:
                    aa = self.codon_to_aa[triplet]
                    self.codon_counts[aa][triplet] += 1
                    self.total_codon_count += 1
                else:
                    # Should theoretically not reach here (pure ATCG should not have non-standard codons)
                    self.skipped_triplets += 1
    
    def calculate_codon_frequencies(self):
        """
        Calculate codon frequencies
        """
        # Calculate totals for each amino acid
        for aa in self.codon_counts:
            self.aa_totals[aa] = sum(self.codon_counts[aa].values())
        
        # Calculate frequency of each codon within synonymous codon group
        for aa in self.codon_counts:
            total = self.aa_totals[aa]
            if total > 0:
                for codon in self.codon_counts[aa]:
                    self.codon_freq[aa][codon] = (self.codon_counts[aa][codon] / total) * 100
    
    def process_sequences(self):
        """
        Process all sequences
        """
        print("Parsing FASTA file...")
        records = self.parse_fasta()
        
        if not records:
            print("No valid sequence records found")
            return
        
        self.total_records = len(records)
        print(f"Found {self.total_records} sequence records")
        
        print("Validating sequence compliance...")
        for header, sequence in records:
            is_valid, violations, warnings, start_codon, stop_codon = self.validate_cds(header, sequence)
            
            if is_valid:
                self.valid_records += 1
                self.valid_lengths.append(len(sequence))
                
                # Further filter compliant sequences: only pure ATCG sequences participate in codon statistics
                if self.is_pure_atcg(sequence):
                    self.analyze_codon_usage(sequence)
                else:
                    # Compliant but non-ATCG sequences are excluded
                    self.excluded_atcg_bad_genes += 1
            else:
                self.invalid_records.append({
                    'gene': header,
                    'length': len(sequence),
                    'start_codon': start_codon,
                    'stop_codon': stop_codon,
                    'violations': '; '.join(violations),
                    'warnings': '; '.join(warnings) if warnings else ""
                })
        
        # Set participating genes count and calculate frequencies
        self.participating_genes = self.valid_records - self.excluded_atcg_bad_genes
        self.calculate_codon_frequencies()
    
    def create_length_histogram(self):
        """
        Create histogram of valid CDS length distribution
        """
        if not self.valid_lengths:
            return ""
        
        # Set Chinese font
        plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
        plt.rcParams['axes.unicode_minus'] = False
        
        plt.figure(figsize=(10, 6))
        plt.hist(self.valid_lengths, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.title('Valid CDS Length Distribution', fontsize=16, fontweight='bold')
        plt.xlabel('CDS Length (bp)', fontsize=12)
        plt.ylabel('Frequency', fontsize=12)
        plt.grid(True, alpha=0.3)
        
        # Add statistical information
        mean_length = sum(self.valid_lengths) / len(self.valid_lengths)
        plt.axvline(mean_length, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_length:.0f} bp')
        plt.legend()
        
        # Save image as base64 string
        buffer = BytesIO()
        plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
        buffer.seek(0)
        image_png = buffer.getvalue()
        buffer.close()
        plt.close()
        
        graphic = base64.b64encode(image_png).decode('utf-8')
        return f'<img src="data:image/png;base64,{graphic}" style="max-width: 100%; height: auto;">'
    
    def generate_html_report(self, output_file="cds_codon_report.html"):
        """
        Generate HTML report
        """
        print("Generating HTML report...")
        
        invalid_count = len(self.invalid_records)
        if self.total_records == 0:
            invalid_percentage = 0.00
            histogram_html = ""
        else:
            invalid_percentage = (invalid_count / self.total_records * 100)
            # Create length histogram
            histogram_html = self.create_length_histogram()
        # ---------- Build data objects for frontend injection ----------
        # 1) CODON_TO_AA
        codon_to_aa_obj = dict(self.codon_to_aa)

        # 2) AA_CODON_FREQ: percentage 0-100, two decimals
        aa_codon_freq_obj = {}
        for aa, codon_map in self.codon_freq.items():
            aa_codon_freq_obj[aa] = {}
            for codon, freq in codon_map.items():
                aa_codon_freq_obj[aa][codon] = round(float(freq), 2)

        # 3) TOP_CODON_PER_AA: top frequency/count within synonymous group; ties by alphabetic order
        top_codon_per_aa_obj = {}
        for aa, codon_map in self.codon_counts.items():
            if not codon_map:
                continue
            # Find maximum count
            max_count = max(codon_map.values())
            candidates = [c for c, cnt in codon_map.items() if cnt == max_count]
            top_codon_per_aa_obj[aa] = sorted(candidates)[0]

        # 4) AA_ONE_LETTER: parse single-letter code from "Ala (A)"; STOP -> '*'
        aa_one_letter_obj = {}
        for aa in set(list(self.codon_counts.keys()) + list(self.aa_totals.keys()) + list(self.codon_freq.keys()) + list(self.codon_to_aa.values())):
            if aa == 'STOP':
                aa_one_letter_obj[aa] = '*'
            else:
                # Extract single-letter code from parentheses
                if '(' in aa and ')' in aa:
                    inside = aa[aa.find('(')+1:aa.find(')')].strip()
                    aa_one_letter_obj[aa] = inside if inside else ''
                else:
                    aa_one_letter_obj[aa] = ''

        # 5) CODON_USAGE_ROWS: table rows in original order
        codon_usage_rows_obj = []
        for aa in sorted(self.aa_totals.keys(), key=lambda a: (a == 'STOP', a)):
            codon_list = [(c, self.codon_counts[aa][c], self.codon_freq[aa][c]) for c in self.codon_counts[aa]]
            codon_list.sort(key=lambda x: (-x[1], x[0]))
            for codon, count, freq in codon_list:
                codon_usage_rows_obj.append({
                    'aa': aa,
                    'codon': codon,
                    'count': count,
                    'freq': round(float(freq), 2)
                })

        # ---------- Number formatting helpers ----------
        def fmt_int(n):
            try:
                return f"{int(n):,}"
            except Exception:
                return str(n)

        def fmt_pct(p):
            try:
                return f"{float(p):.2f}%"
            except Exception:
                return str(p)

        invalid_preview_limit = 3

        def render_invalid_row(r):
            return (
                f"            <tr><td>{html.escape(r['gene'])}</td><td>{fmt_int(r['length'])}</td>"
                f"<td>{html.escape(r['start_codon'])}</td><td>{html.escape(r['stop_codon'])}</td>"
                f"<td class='violations'>{html.escape(r['violations'])}</td>"
                f"<td class='warnings'>{html.escape(r['warnings'])}</td></tr>\n"
            )

        # ---------- Styles and base scripts (non f-string) ----------
        css_styles = """
        :root {
            --bg-page: #f3f6fb;
            --bg-panel: #ffffff;
            --bg-subtle: #f6f8fc;
            --bg-accent-soft: #eef5ff;
            --text-primary: #1f2d3d;
            --text-secondary: #5f6c7a;
            --border-color: #d8e0ea;
            --primary: #2a6fdb;
            --primary-dark: #1f5abc;
            --success: #1f9d67;
            --danger: #cc3d3d;
            --warning: #c27c0e;
            --shadow-soft: 0 10px 24px rgba(31, 45, 61, 0.08);
            --radius-md: 10px;
            --radius-sm: 8px;
        }
        body {
            font-family: "Segoe UI", "Helvetica Neue", Arial, sans-serif;
            margin: 20px;
            background-color: var(--bg-page);
            color: var(--text-primary);
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background-color: var(--bg-panel);
            padding: 24px;
            border-radius: 14px;
            box-shadow: var(--shadow-soft);
        }
        h1 {
            color: var(--text-primary);
            text-align: center;
            border-bottom: 3px solid var(--primary);
            padding-bottom: 12px;
            letter-spacing: 0.2px;
        }
        h2 {
            color: var(--text-primary);
            border-left: 4px solid var(--primary);
            padding-left: 10px;
            margin-top: 32px;
        }
        .summary {
            background-color: var(--bg-subtle);
            border: 1px solid var(--border-color);
            padding: 16px;
            border-radius: var(--radius-md);
            margin: 20px 0;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
            gap: 12px;
        }
        .summary-card {
            background: var(--bg-panel);
            border: 1px solid var(--border-color);
            border-radius: var(--radius-sm);
            padding: 12px;
            box-shadow: 0 2px 8px rgba(31, 45, 61, 0.05);
        }
        .summary-label {
            color: var(--text-secondary);
            font-size: 13px;
            margin-bottom: 6px;
        }
        .summary-value {
            color: var(--text-primary);
            font-size: 22px;
            font-weight: 700;
            line-height: 1.2;
        }
        .summary-value.valid { color: var(--success); }
        .summary-value.invalid { color: var(--danger); }
        .summary-warning {
            color: var(--danger);
            font-weight: 700;
            margin-top: 10px;
        }
        .summary-item {
            display: inline-block;
            margin: 8px 16px 8px 0;
            color: var(--text-primary);
            font-weight: 600;
        }
        .valid { color: var(--success); }
        .invalid { color: var(--danger); }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 14px;
            table-layout: fixed;
        }
        th, td {
            border: 1px solid var(--border-color);
            padding: 10px 8px;
            text-align: left;
            word-wrap: break-word;
            word-break: break-all;
            overflow-wrap: break-word;
            line-height: 1.4;
        }
        .gene-column { width: 35%; max-width: 35%; }
        .length-column { width: 10%; max-width: 10%; }
        .codon-column { width: 12%; max-width: 12%; }
        .violations-column { width: 22%; max-width: 22%; }
        .warnings-column { width: 15%; max-width: 15%; }
        .non-atcg-column { width: 57%; max-width: 57%; }
        .codon-aa-column { width: 25%; max-width: 25%; }
        .codon-freq-column { width: 18%; max-width: 18%; }
        th {
            background-color: var(--primary);
            color: #ffffff;
            font-weight: bold;
        }
        tr:nth-child(even) { background-color: #f9fbff; }
        tr:hover { background-color: var(--bg-accent-soft); }
        .violations { color: var(--danger); font-weight: bold; }
        .warnings { color: var(--warning); font-style: italic; }
        .meta {
            background-color: #f8f9fa;
            padding: 15px;
            border-radius: var(--radius-sm);
            border: 1px solid var(--border-color);
            margin-top: 30px;
            font-size: 12px;
            color: var(--text-secondary);
        }
        .violation-stats {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 12px;
            margin: 20px 0;
        }
        .violation-item {
            text-align: left;
            background-color: var(--bg-panel);
            padding: 14px;
            border-radius: var(--radius-sm);
            border: 1px solid var(--border-color);
            box-shadow: 0 2px 8px rgba(31, 45, 61, 0.05);
        }
        .violation-number { font-size: 24px; font-weight: bold; color: var(--danger); }
        .violation-label {
            margin-top: 6px;
            color: var(--text-secondary);
            font-size: 13px;
        }
        .violation-rate {
            margin-top: 2px;
            color: var(--text-primary);
            font-weight: 600;
            font-size: 13px;
        }
        .histogram-container { text-align: center; margin: 20px 0; }
        .invalid-records-preview-note {
            margin: 8px 0 12px 0;
            color: var(--text-secondary);
            font-size: 14px;
        }
        .invalid-records-details {
            margin-top: 8px;
            border: 1px solid var(--border-color);
            border-radius: var(--radius-sm);
            background: var(--bg-subtle);
            padding: 8px 10px;
        }
        .invalid-records-details summary {
            cursor: pointer;
            font-weight: bold;
            color: var(--primary-dark);
            margin: 2px 0 8px 0;
            list-style: none;
        }
        .invalid-records-details summary::-webkit-details-marker {
            display: none;
        }
        .invalid-records-details summary::before {
            content: "▸ ";
            color: var(--primary-dark);
            font-weight: 700;
        }
        .invalid-records-details[open] summary::before {
            content: "▾ ";
        }
        /* Upload & Optimize base styles */
        .upload-card {
            background: linear-gradient(180deg, #f9fbff 0%, #f4f8ff 100%);
            border: 1px solid var(--border-color);
            padding: 18px;
            border-radius: var(--radius-md);
            margin: 20px 0;
            box-shadow: 0 6px 16px rgba(31, 45, 61, 0.06);
        }
        .upload-card p { margin: 0 0 12px 0; color: var(--text-secondary); line-height: 1.45; }
        .controls { display: flex; gap: 10px; flex-wrap: wrap; align-items: center; }
        .btn {
            background: var(--primary);
            color: #fff;
            border: none;
            padding: 9px 14px;
            border-radius: 4px;
            cursor: pointer;
            transition: background-color .2s ease;
            font-weight: 600;
        }
        .btn:disabled { cursor: not-allowed; opacity: 0.6; }
        .btn:hover:not(:disabled) { background: var(--primary-dark); }
        .btn.btn-secondary {
            background: #ffffff;
            color: var(--text-secondary);
            border: 1px solid var(--border-color);
        }
        .btn.btn-secondary:hover:not(:disabled) {
            background: #f2f5fa;
            color: var(--text-primary);
        }
        #btnOptimize {
            background: linear-gradient(180deg, #2f78ec 0%, #2668d0 100%);
            box-shadow: 0 4px 10px rgba(38, 104, 208, 0.3);
        }
        #btnOptimize:hover:not(:disabled) {
            background: linear-gradient(180deg, #286adb 0%, #1f5abc 100%);
        }
        .results {
            background: #ecf0f1;
            border: 1px dashed #bdc3c7;
            padding: 10px;
            border-radius: 6px;
            margin-top: 12px;
            display: none;
        }
        .replaced { color: #c0392b; font-weight: bold; }
        .toast { display:none; padding:10px; border-radius:6px; margin-top:8px; }
        .toast[aria-live] { outline: none; }
        .toast-success { background:#e8f6ee; border:1px solid #27ae60; color:#1e824c; }
        .toast-error { background:#fdecea; border:1px solid #e74c3c; color:#c0392b; }
        .toast-info { background:#eef5ff; border:1px solid #2980b9; color:#2c3e50; }
        .result-success { background: #e8f6ee; border: 1px solid #27ae60; color: #1e824c; }
        .result-error { background: #fdecea; border: 1px solid #e74c3c; color: #c0392b; }
        .mono-box {
            background: #ffffff;
            border: 1px solid #bdc3c7;
            padding: 10px;
            border-radius: 6px;
            font-family: "Courier New", Courier, monospace;
            overflow-x: auto;
            white-space: pre-wrap;
        }
        /* Codon Usage Table Enhancements */
        .codon-usage-filters {
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            padding: 15px;
            margin: 20px 0;
        }
        .filter-controls {
            display: flex;
            align-items: center;
            gap: 15px;
            flex-wrap: wrap;
            margin-bottom: 10px;
        }
        .filter-label {
            font-weight: bold;
            color: #495057;
            white-space: nowrap;
        }
        .filter-select, .filter-input {
            padding: 6px 10px;
            border: 1px solid #ced4da;
            border-radius: 4px;
            font-size: 14px;
        }
        .filter-select {
            min-width: 120px;
        }
        .filter-input {
            min-width: 200px;
        }
        .filter-btn {
            padding: 6px 12px;
            font-size: 14px;
        }
        .row-count {
            font-weight: bold;
            color: #495057;
            margin-left: auto;
        }
        .export-controls {
            display: flex;
            align-items: center;
            gap: 15px;
            flex-wrap: wrap;
            border-top: 1px solid #dee2e6;
            padding-top: 10px;
        }
        .export-checkbox-label {
            display: flex;
            align-items: center;
            gap: 5px;
            font-size: 14px;
            color: #495057;
        }
        .export-checkbox {
            margin: 0;
        }
        .export-btn {
            padding: 6px 12px;
            font-size: 14px;
        }
        .export-btn:disabled {
            opacity: 0.6;
            cursor: not-allowed;
        }
        .export-notes {
            margin-top: 10px;
            color: #6c757d;
            font-size: 12px;
            line-height: 1.4;
        }
        .codon-usage-table-wrapper {
            max-height: 520px;
            overflow: auto;
            border: 1px solid #dee2e6;
            border-radius: 6px;
            margin: 20px 0;
        }
        #codonUsageTable {
            margin: 0;
            border-collapse: collapse;
        }
        #codonUsageTable thead th {
            position: sticky;
            top: var(--sticky-offset, 64px);
            z-index: 2;
            background-color: #3498db;
            color: white;
            font-weight: bold;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        #codonUsageTable tbody tr.group-start {
            border-top: 2px solid #2c3e50;
        }
        #codonUsageTable tbody tr.group-start td {
            border-top: 2px solid #2c3e50;
        }
        .report-nav {
            position: sticky;
            top: 8px;
            z-index: 30;
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
            background: rgba(255, 255, 255, 0.92);
            backdrop-filter: blur(6px);
            border: 1px solid var(--border-color);
            border-radius: var(--radius-sm);
            padding: 8px;
            margin: 12px 0 16px 0;
            box-shadow: 0 6px 16px rgba(31, 45, 61, 0.08);
        }
        .nav-link {
            display: inline-block;
            text-decoration: none;
            color: var(--primary-dark);
            background: #f4f8ff;
            border: 1px solid #dce7fb;
            border-radius: 999px;
            padding: 6px 10px;
            font-size: 13px;
            font-weight: 600;
        }
        .nav-link:hover {
            background: #e9f1ff;
            border-color: #c9dcff;
        }
        .nav-link:focus-visible {
            outline: 2px solid var(--primary);
            outline-offset: 2px;
        }
        .report-section {
            background: var(--bg-panel);
            border: 1px solid var(--border-color);
            border-radius: var(--radius-md);
            padding: 16px;
            margin-top: 16px;
            scroll-margin-top: 88px;
            box-shadow: 0 4px 14px rgba(31, 45, 61, 0.05);
        }
        .report-section h2 {
            margin-top: 0;
        }
        .section-note {
            color: var(--text-secondary);
            margin: 4px 0 10px 0;
            font-size: 13px;
        }
        .back-to-top {
            position: fixed;
            right: 16px;
            bottom: 16px;
            z-index: 40;
            border: none;
            border-radius: 999px;
            background: var(--primary-dark);
            color: #fff;
            padding: 10px 14px;
            font-size: 13px;
            font-weight: 700;
            cursor: pointer;
            box-shadow: 0 8px 20px rgba(31, 90, 188, 0.35);
            opacity: 0;
            transform: translateY(8px);
            pointer-events: none;
            transition: opacity .2s ease, transform .2s ease;
        }
        .back-to-top.is-visible {
            opacity: 1;
            transform: translateY(0);
            pointer-events: auto;
        }
        .back-to-top:focus-visible {
            outline: 2px solid #ffffff;
            outline-offset: 2px;
        }
        @media (max-width: 768px) {
            .container { padding: 10px; margin: 10px; }
            .report-nav { top: 4px; padding: 6px; }
            .nav-link { font-size: 12px; padding: 6px 8px; }
            .report-section { padding: 12px; }
            .back-to-top { right: 12px; bottom: 12px; }
            .summary-grid { grid-template-columns: 1fr 1fr; }
            .summary-value { font-size: 20px; }
            .violation-stats { grid-template-columns: 1fr; }
            table { font-size: 12px; }
            th, td { padding: 6px 4px; }
            .gene-column { width: 40%; max-width: 40%; }
            .length-column { width: 12%; max-width: 12%; }
            .codon-column { width: 12%; max-width: 12%; }
            .violations-column { width: 24%; max-width: 24%; }
            .warnings-column { width: 12%; max-width: 12%; }
            .non-atcg-column { width: 48%; max-width: 48%; }
            .codon-aa-column { width: 30%; max-width: 30%; }
            .codon-freq-column { width: 18%; max-width: 18%; }
            .controls { flex-direction: column; align-items: stretch; }
            .controls .btn { width: 100%; }
        }
        """

        js_scripts = r"""
        (function(){
            // Codon Usage Table Management
            var codonUsageRows = window.CODON_USAGE_ROWS || [];
            var filteredRows = codonUsageRows;
            var aaFilter = document.getElementById('aaFilter');
            var searchFilter = document.getElementById('searchFilter');
            var clearFilters = document.getElementById('clearFilters');
            var rowCount = document.getElementById('rowCount');
            var exportFullDataset = document.getElementById('exportFullDataset');
            var downloadXLSX = document.getElementById('downloadCodonUsageXLSX');
            var downloadCSV = document.getElementById('downloadCodonUsageCSV');
            var tableBody = document.getElementById('codonUsageTableBody');
            var searchTimeout = null;
            var reportNav = document.getElementById('reportNav');
            var backToTopBtn = document.getElementById('backToTopBtn');

            function updateStickyOffset() {
                var navHeight = reportNav ? reportNav.offsetHeight : 0;
                document.documentElement.style.setProperty('--sticky-offset', (navHeight + 14) + 'px');
            }

            function handleBackToTopVisibility() {
                if (!backToTopBtn) return;
                if (window.scrollY > 260) {
                    backToTopBtn.classList.add('is-visible');
                } else {
                    backToTopBtn.classList.remove('is-visible');
                }
            }

            function initPageNavigation() {
                updateStickyOffset();
                handleBackToTopVisibility();
                window.addEventListener('resize', updateStickyOffset);
                window.addEventListener('scroll', handleBackToTopVisibility, { passive: true });
                if (backToTopBtn) {
                    backToTopBtn.addEventListener('click', function() {
                        window.scrollTo({ top: 0, behavior: 'smooth' });
                    });
                }
            }

            // Initialize codon usage table
            function initCodonUsageTable() {
                if (!codonUsageRows.length || !tableBody) return;
                
                // Populate amino acid filter
                if (aaFilter) {
                    var uniqueAAs = [...new Set(codonUsageRows.map(row => row.aa))].sort();
                    uniqueAAs.forEach(aa => {
                        var option = document.createElement('option');
                        option.value = aa;
                        option.textContent = aa;
                        aaFilter.appendChild(option);
                    });
                }
                
                // Render initial table
                renderCodonUsageTable();
                
                // Add event listeners
                if (aaFilter) aaFilter.addEventListener('change', applyFilters);
                if (searchFilter) {
                    searchFilter.addEventListener('input', function() {
                        clearTimeout(searchTimeout);
                        searchTimeout = setTimeout(applyFilters, 200);
                    });
                }
                if (clearFilters) clearFilters.addEventListener('click', clearAllFilters);
                if (downloadXLSX) downloadXLSX.addEventListener('click', downloadXLSXReport);
                if (downloadCSV) downloadCSV.addEventListener('click', downloadCSVReport);
            }

            function applyFilters() {
                var aaValue = aaFilter ? aaFilter.value : '';
                var searchValue = searchFilter ? searchFilter.value.toLowerCase() : '';
                
                filteredRows = codonUsageRows.filter(function(row) {
                    var aaMatch = !aaValue || row.aa === aaValue;
                    var searchMatch = !searchValue || 
                        row.aa.toLowerCase().includes(searchValue) || 
                        row.codon.toLowerCase().includes(searchValue);
                    return aaMatch && searchMatch;
                });
                
                renderCodonUsageTable();
                updateRowCount();
                updateExportButtons();
            }

            function clearAllFilters() {
                if (aaFilter) aaFilter.value = '';
                if (searchFilter) searchFilter.value = '';
                filteredRows = codonUsageRows;
                renderCodonUsageTable();
                updateRowCount();
                updateExportButtons();
            }

            function renderCodonUsageTable() {
                if (!tableBody) return;
                
                var html = '';
                var lastAA = '';
                
                filteredRows.forEach(function(row, index) {
                    var isGroupStart = row.aa !== lastAA;
                    var groupClass = isGroupStart ? ' group-start' : '';
                    
                    html += '<tr class="' + groupClass + '">';
                    html += '<td>' + escapeHtml(row.aa) + '</td>';
                    html += '<td>' + escapeHtml(row.codon) + '</td>';
                    html += '<td>' + formatNumber(row.count) + '</td>';
                    html += '<td>' + row.freq.toFixed(2) + '%</td>';
                    html += '</tr>';
                    
                    lastAA = row.aa;
                });
                
                tableBody.innerHTML = html;
            }

            function updateRowCount() {
                if (rowCount) {
                    rowCount.textContent = 'Rows: ' + filteredRows.length + ' / ' + codonUsageRows.length;
                }
            }

            function updateExportButtons() {
                var hasData = filteredRows.length > 0;
                var disabled = !hasData;
                
                if (downloadXLSX) {
                    downloadXLSX.disabled = disabled;
                    downloadXLSX.title = disabled ? 'No codon usage data to export (filtered out)' : 'Export filtered codon usage data as XLSX';
                }
                if (downloadCSV) {
                    downloadCSV.disabled = disabled;
                    downloadCSV.title = disabled ? 'No codon usage data to export (filtered out)' : 'Export filtered codon usage data as CSV';
                }
            }

            function downloadXLSXReport() {
                var data = exportFullDataset && exportFullDataset.checked ? codonUsageRows : filteredRows;
                if (!data.length) return;
                
                generateXLSXReport(data);
            }

            function downloadCSVReport() {
                var data = exportFullDataset && exportFullDataset.checked ? codonUsageRows : filteredRows;
                if (!data.length) return;
                
                var csv = generateCSVReport(data);
                downloadFile(csv, 'Codon_usage.csv', 'text/csv;charset=utf-8;');
            }

            function generateXLSXReport(data) {
                // Create a new workbook and worksheet
                var workbook = new ExcelJS.Workbook();
                var worksheet = workbook.addWorksheet('CodonUsage');
                
                // Set column headers
                worksheet.columns = [
                    { header: 'Amino acid', key: 'aa', width: 20 },
                    { header: 'Codon', key: 'codon', width: 10 },
                    { header: 'Count', key: 'count', width: 16 },
                    { header: 'Frequency within AA (%)', key: 'freq', width: 24 }
                ];
                
                // Add data rows
                var lastAA = '';
                data.forEach(function(row, index) {
                    var isGroupStart = row.aa !== lastAA;
                    
                    // Convert frequency from percentage string to decimal if needed
                    var freqValue = row.freq;
                    if (typeof freqValue === 'string' && freqValue.includes('%')) {
                        freqValue = parseFloat(freqValue.replace('%', '')) / 100;
                    } else if (typeof freqValue === 'number' && freqValue > 1) {
                        freqValue = freqValue / 100;
                    }
                    
                    var rowData = {
                        aa: row.aa,
                        codon: row.codon,
                        count: row.count,
                        freq: freqValue
                    };
                    
                    var excelRow = worksheet.addRow(rowData);
                    
                    // Add thick top border for group start rows
                    if (isGroupStart && index > 0) {
                        excelRow.getCell(1).border = {
                            top: { style: 'thick', color: { argb: 'FF000000' } }
                        };
                        excelRow.getCell(2).border = {
                            top: { style: 'thick', color: { argb: 'FF000000' } }
                        };
                        excelRow.getCell(3).border = {
                            top: { style: 'thick', color: { argb: 'FF000000' } }
                        };
                        excelRow.getCell(4).border = {
                            top: { style: 'thick', color: { argb: 'FF000000' } }
                        };
                    }
                    
                    lastAA = row.aa;
                });
                
                // Set number formats
                worksheet.getColumn('count').numFmt = '#,##0';
                worksheet.getColumn('freq').numFmt = '0.00%';
                
                // Freeze the header row
                worksheet.views = [
                    { state: 'frozen', ySplit: 1 }
                ];
                
                // Enable AutoFilter
                var lastRow = worksheet.rowCount;
                worksheet.autoFilter = 'A1:D' + lastRow;
                
                // Generate and download the file
                workbook.xlsx.writeBuffer().then(function(buffer) {
                    var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
                    saveAs(blob, 'Codon_usage.xlsx');
                }).catch(function(error) {
                    console.error('Error generating XLSX:', error);
                    alert('Error generating XLSX file. Please try again.');
                });
            }

            function generateCSVReport(data) {
                var csv = 'Amino acid,Codon,Count,Frequency within AA (%)\n';
                data.forEach(function(row) {
                    // Normalize Count only for CSV export: if string, strip thousands separators
                    var countForCsv;
                    if (typeof row.count === 'string') {
                        countForCsv = row.count.replace(/,/g, '');
                    } else {
                        countForCsv = String(row.count);
                    }
                    // Keep current frequency export format as-is (two decimals, no extra symbols here)
                    var freqValue = row.freq.toFixed(2);
                    csv += csvCell(row.aa) + ',' + csvCell(row.codon) + ',' + csvCell(countForCsv) + ',' + csvCell(freqValue) + '\r\n';
                });
                return csv;
            }

            function downloadFile(content, filename, mimeType) {
                // Add BOM for better Excel compatibility when downloading CSV
                var finalContent = content;
                if (mimeType === 'text/csv' || mimeType === 'text/csv;charset=utf-8;') {
                    finalContent = '\uFEFF' + content;
                }
                
                var blob = new Blob([finalContent], {type: mimeType});
                var url = URL.createObjectURL(blob);
                var a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                URL.revokeObjectURL(url);
            }

            function escapeHtml(text) {
                var div = document.createElement('div');
                div.textContent = text;
                return div.innerHTML;
            }

            function csvCell(value) {
                // Convert to string
                var s = value.toString();
                
                // RFC4180: If contains comma, quote, newline, or carriage return, wrap in quotes and escape internal quotes
                if (s.includes(',') || s.includes('"') || s.includes('\n') || s.includes('\r')) {
                    return '"' + s.replace(/"/g, '""') + '"';
                }
                
                return s;
            }
            
            function escapeCsv(text) {
                if (text.includes(',') || text.includes('"') || text.includes('\n')) {
                    return '"' + text.replace(/"/g, '""') + '"';
                }
                return text;
            }

            function formatNumber(num) {
                return num.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ',');
            }

            // Initialize codon usage table when DOM is ready
            if (document.readyState === 'loading') {
                document.addEventListener('DOMContentLoaded', function() {
                    initCodonUsageTable();
                    initPageNavigation();
                });
            } else {
                initCodonUsageTable();
                initPageNavigation();
            }

            // Original Upload & Optimize functionality
            var fileInput = document.getElementById('faUpload');
            var btnOptimize = document.getElementById('btnOptimize');
            var btnClear = document.getElementById('btnClear');
            var uploadResults = document.getElementById('uploadResults');
            var optimizeResults = document.getElementById('optimizeResults');
            var statusToast = document.getElementById('statusToast');

            var state = { compliant: false, cdsName: '', seq: '', filename: '' };

            function showToast(text, type){
                if (!statusToast) return;
                statusToast.textContent = text || '';
                statusToast.className = 'toast ' + (type || 'toast-info');
                statusToast.style.display = text ? 'block' : 'none';
            }

            function resetUI(){
                state = { compliant: false, cdsName: '', seq: '', filename: '' };
                if (fileInput) fileInput.value = '';
                if (uploadResults) { uploadResults.style.display = 'none'; uploadResults.className = 'results'; uploadResults.innerHTML = ''; }
                if (optimizeResults) { optimizeResults.style.display = 'none'; optimizeResults.className = 'results'; optimizeResults.innerHTML = ''; }
                if (btnOptimize) btnOptimize.disabled = true;
                showToast('', 'toast-info');
            }

            function hasFaExtension(name){
                return typeof name === 'string' && name.toLowerCase().endsWith('.fa');
            }

            function parseFastaSingle(text){
                var lines = (text || '').split(/\r?\n/);
                var headers = [];
                for (var i=0;i<lines.length;i++){
                    if (lines[i].trim().startsWith('>')) headers.push(i);
                }
                if (headers.length !== 1) {
                    return { error: 'FASTA must contain exactly one record (single header line starting with ">").' };
                }
                var headerIdx = headers[0];
                var headerLine = lines[headerIdx].slice(1).trim();
                var seqParts = [];
                for (var j=headerIdx+1;j<lines.length;j++){
                    var ln = lines[j].trim();
                    if (!ln) continue;
                    if (ln.startsWith('>')){
                        return { error: 'Multiple FASTA records detected. Please upload only one record.' };
                    }
                    seqParts.push(ln);
                }
                var seq = seqParts.join('').replace(/\s+/g,'').toUpperCase();
                if (!seq) return { error: 'Sequence is empty.' };
                return { header: headerLine, sequence: seq };
            }

            function validateSeq(seq){
                var reasons = [];
                var nonATCG = [];
                var validSet = { 'A':true, 'T':true, 'C':true, 'G':true };
                var setSeen = {};
                for (var i=0;i<seq.length;i++){
                    var ch = seq[i];
                    if (!validSet[ch]) setSeen[ch] = true;
                }
                nonATCG = Object.keys(setSeen).sort();
                if (nonATCG.length>0){
                    reasons.push('Only ATCG check failed: illegal chars [' + nonATCG.join(', ') + ']');
                }

                var lenOk = (seq.length % 3) === 0;
                if (!lenOk) reasons.push('Length not multiple of 3.');

                var startOk = false;
                if (seq.length >= 3){
                    startOk = (seq.slice(0,3) === 'ATG');
                    if (!startOk) reasons.push('Start codon must be ATG.');
                } else {
                    reasons.push('Start codon must be ATG.');
                }
                var missingStop = false;
                if (lenOk && seq.length >= 3){
                    var stop = seq.slice(-3);
                    var stopOk = (stop === 'TGA' || stop === 'TAA' || stop === 'TAG');
                    if (!stopOk) missingStop = true;
                }

                return { ok: reasons.length===0, reasons: reasons, nonATCG: nonATCG, missingStop: missingStop };
            }

            function translateProtein(seq){
                var codonToAA = window.CODON_TO_AA || {};
                var aaOne = window.AA_ONE_LETTER || {};
                var protein = '';
                for (var i=0;i<seq.length; i+=3){
                    if (i+3>seq.length) break;
                    var codon = seq.slice(i,i+3);
                    var aaName = codonToAA[codon];
                    if (!aaName){ protein += 'X'; continue; }
                    var letter = aaOne[aaName];
                    protein += (letter && typeof letter==='string' ? letter : 'X');
                }
                return protein;
            }

            function wrap60(s){
                var out = [];
                for (var i=0;i<s.length;i+=60){ out.push(s.slice(i, i+60)); }
                return out.join('\n');
            }

            function renderCompliant(filename, header, seq, opts){
                var protein = translateProtein(seq);
                var fastaHeader = '>Translated_protein_of_' + header;
                var proteinFasta = fastaHeader + '\n' + wrap60(protein);
                var notes = '';
                if (opts && opts.missingStop){
                    notes = '<div class="results" style="display:block;background:#fffbea;border:1px solid #f1c40f;color:#8a6d3b;margin-top:10px;">Note: no STOP codon detected (allowed).</div>';
                }
                var html = ''+
                    '<div><strong>File:</strong> ' + (filename || '') + '</div>'+
                    '<div><strong>CDS name:</strong> ' + header + '</div>'+
                    '<div><strong>CDS length:</strong> ' + seq.length.toLocaleString() + ' bp</div>'+
                    '<div class="results result-success" style="display:block;margin-top:10px;">✅ This CDS is compliant.</div>'+
                    notes +
                    '<div style="margin-top:10px;"><strong>Translated protein (FASTA):</strong></div>'+
                    '<div class="mono-box">' + proteinFasta.replace(/&/g,'&amp;').replace(/</g,'&lt;') + '</div>';
                uploadResults.innerHTML = html;
                uploadResults.style.display = 'block';
                uploadResults.className = 'results';
                state.compliant = true;
                state.protein = protein;
                btnOptimize.disabled = false;
            }

            function renderErrors(reasons){
                var list = '<ul style="margin:0;padding-left:18px;">' + reasons.map(function(r){ return '<li>'+r.replace(/&/g,'&amp;').replace(/</g,'&lt;')+'</li>'; }).join('') + '</ul>';
                uploadResults.innerHTML = '<div class="results result-error" style="display:block;">'+ list + '</div>';
                uploadResults.style.display = 'block';
                uploadResults.className = 'results';
                state.compliant = false;
                btnOptimize.disabled = true;
            }

            function handleFileChange(){
                optimizeResults.style.display = 'none';
                optimizeResults.innerHTML = '';
                optimizeResults.className = 'results';
                var file = (fileInput && fileInput.files && fileInput.files[0]) ? fileInput.files[0] : null;
                if (!file){ resetUI(); return; }
                if (!hasFaExtension(file.name)){
                    resetUI();
                    showToast('Only .fa extension is accepted.', 'toast-error');
                    return;
                }
                state.filename = file.name;
                var reader = new FileReader();
                reader.onload = function(ev){
                    var text = ev.target.result || '';
                    var parsed = parseFastaSingle(text);
                    if (parsed.error){
                        showToast(parsed.error, 'toast-error');
                        renderErrors([parsed.error]);
                        return;
                    }
                    state.cdsName = parsed.header;
                    state.seq = parsed.sequence;
                    var v = validateSeq(parsed.sequence);
                    if (v.ok){
                        showToast('File loaded and validated.', 'toast-success');
                        renderCompliant(state.filename, state.cdsName, state.seq, { missingStop: v.missingStop });
                    } else {
                        showToast('Validation failed. Please review the errors.', 'toast-error');
                        renderErrors(v.reasons);
                    }
                };
                reader.onerror = function(){
                    showToast('Failed to read the file.', 'toast-error');
                    renderErrors(['Failed to read the file.']);
                };
                reader.readAsText(file);
            }

            if (fileInput) fileInput.addEventListener('change', handleFileChange);
            if (btnClear) btnClear.addEventListener('click', function(){ resetUI(); });

            function optimizeCodons(seq){
                var codonToAA = window.CODON_TO_AA || {};
                var topPer = window.TOP_CODON_PER_AA || {};
                var COL_W = 4; // fixed column width, ensure vertical alignment (accommodate STOP)
                function padRight(s, w){ s = String(s); if (s.length >= w) return s; return s + Array(w - s.length + 1).join(' '); }
                function padLeft(s, w){ s = String(s); if (s.length >= w) return s; return Array(w - s.length + 1).join(' ') + s; }
                var parts = [];
                var codonView = [];
                var aaView = [];
                var idxView = [];
                var replaced = 0;
                var total = 0;
                for (var i=0;i<seq.length;i+=3){
                    if (i+3>seq.length) break;
                    var codon = seq.slice(i,i+3);
                    total += 1;
                    var aa = codonToAA[codon];
                    var top = aa ? topPer[aa] : undefined;
                    var use = (aa && top) ? top : codon;
                    var codonCell = codon;
                    if (aa && top && top !== codon){
                        replaced += 1;
                        codonCell = '<span class="replaced">'+codon+'</span>';
                    }
                    // pad with spaces to align columns visually
                    codonView.push(codonCell + padRight('', COL_W - 3));
                    parts.push(use);

                    // AA three-letter label
                    var aaLabel = '';
                    if (aa) {
                        if (aa === 'STOP') { aaLabel = 'STOP'; }
                        else {
                            var sp = aa.indexOf(' ');
                            aaLabel = sp > 0 ? aa.substring(0, sp) : aa;
                        }
                    } else { aaLabel = '?'; }
                    aaView.push(padRight(aaLabel, COL_W));
                    idxView.push(padRight(String(total), COL_W));
                }
                var viewHtml = [
                    codonView.join('·'),
                    aaView.join('·'),
                    idxView.join('·')
                ].join('\n');
                return { optimized: parts.join(''), replaced: replaced, total: total, viewHtml: viewHtml };
            }

            function copyText(text){
                try {
                    navigator.clipboard.writeText(text);
                } catch(e) {
                    var ta = document.createElement('textarea');
                    ta.value = text; document.body.appendChild(ta); ta.select();
                    try{ document.execCommand('copy'); }catch(_){ }
                    document.body.removeChild(ta);
                }
            }

            function safeFileName(name){
                if (!name) return 'unknown';
                var cleaned = String(name).trim().replace(/[^A-Za-z0-9._-]+/g, '_');
                if (!cleaned) cleaned = 'unknown';
                if (cleaned.length > 100) cleaned = cleaned.slice(0, 100);
                return cleaned;
            }

            function renderOptimize(){
                if (!state.compliant || !state.seq) return;
                var res = optimizeCodons(state.seq);
                var optProt = translateProtein(res.optimized);
                var preserved = (optProt === (state.protein || translateProtein(state.seq)));
                var header = '>Optimized_CDS_of_' + (state.cdsName || 'unknown');
                var fasta = header + '\n' + wrap60(res.optimized);
                var translateCheck = preserved ? '✅ preserved' : '❌ changed';
                var actions = ''+
                    '<div class="controls" style="margin-top:8px;gap:8px;">'
                    + '<button id="btnCopyFasta" class="btn" aria-label="Copy optimized CDS to clipboard">Copy optimized CDS</button>'
                    + (preserved ? '<button id="btnDownloadFasta" class="btn" aria-label="Download optimized CDS as FASTA">Download optimized CDS (.fa)</button>' : '<button id="btnDownloadFasta" class="btn" disabled aria-label="Download disabled due to translation mismatch">Download optimized CDS (.fa)</button>')
                    + '</div>';
                var warn = preserved ? '' : '<div class="results result-error" style="display:block;margin-top:8px;">Translation changed after optimization. Please review.</div>';
                var html = ''+
                    '<div><strong>Replacements:</strong> ' + res.replaced.toLocaleString() + ' / ' + res.total.toLocaleString() + '</div>'+
                    '<div style="margin-top:6px;"><strong>Translation check:</strong> ' + translateCheck + '</div>'+
                    '<div style="margin-top:8px;"><strong>Codon-wise view:</strong></div>'+
                    '<div class="mono-box" style="line-height:1.6; white-space:pre;">' + res.viewHtml + '</div>'+
                    warn +
                    '<div style="margin-top:10px;"><strong>Optimized CDS (FASTA):</strong></div>'+
                    '<div class="mono-box">' + fasta.replace(/&/g,'&amp;').replace(/</g,'&lt;') + '</div>'+
                    actions;
                optimizeResults.innerHTML = html;
                optimizeResults.style.display = 'block';
                optimizeResults.className = 'results';
                var btnCopy = document.getElementById('btnCopyFasta');
                if (btnCopy){ btnCopy.addEventListener('click', function(){ copyText(fasta); showToast('Optimized FASTA has been copied to clipboard.', 'toast-success'); }); }
                var btnDl = document.getElementById('btnDownloadFasta');
                if (btnDl && preserved){
                    btnDl.addEventListener('click', function(){
                        var clean = safeFileName(state.cdsName || 'unknown');
                        var filename = 'optimized_' + clean + '.fa';
                        try {
                            var blob = new Blob([fasta], {type: 'text/plain;charset=utf-8'});
                            var url = URL.createObjectURL(blob);
                            var a = document.createElement('a');
                            a.href = url; a.download = filename;
                            document.body.appendChild(a);
                            a.click();
                            document.body.removeChild(a);
                            setTimeout(function(){ URL.revokeObjectURL(url); }, 0);
                            showToast('Download started: ' + filename, 'toast-success');
                        } catch (e) {
                            showToast('Failed to trigger download.', 'toast-error');
                        }
                    });
                }
            }

            if (btnOptimize) btnOptimize.addEventListener('click', renderOptimize);

            window.__CDS_REPORT_READY__ = true;
        })();
        """

        # ---------- Pre-construct major HTML sections to avoid backslashes inside f-strings ----------
        # Invalid Records
        if self.invalid_records:
            invalid_table_header = """
        <table>
            <tr>
                <th class="gene-column">Gene</th>
                <th class="length-column">Length</th>
                <th class="codon-column">Start codon</th>
                <th class="codon-column">Stop codon</th>
                <th class="violations-column">Violations</th>
                <th class="warnings-column">Warnings</th>
            </tr>
"""
            preview_rows = ''.join(render_invalid_row(r) for r in self.invalid_records[:invalid_preview_limit])

            if len(self.invalid_records) > invalid_preview_limit:
                hidden_rows = ''.join(render_invalid_row(r) for r in self.invalid_records[invalid_preview_limit:])
                hidden_count = len(self.invalid_records) - invalid_preview_limit
                invalid_section_html = (
                    f"<div class='invalid-records-preview-note'>Showing first {invalid_preview_limit} records. Total invalid records: {fmt_int(len(self.invalid_records))}.</div>"
                    + invalid_table_header
                    + preview_rows
                    + "        </table>\n"
                    + f"<details class='invalid-records-details'><summary>Show remaining {fmt_int(hidden_count)} invalid records</summary>"
                    + invalid_table_header
                    + hidden_rows
                    + "        </table></details>\n"
                )
            else:
                invalid_section_html = invalid_table_header + preview_rows + "        </table>\n"
        else:
            invalid_section_html = "<p style=\"color: #27ae60; font-weight: bold;\">🎉 All sequences are compliant!</p>"

        # Non-ATCG Records
        if self.non_atcg_records:
            nonatcg_table_header = """
        <table>
            <tr>
                <th class="gene-column">Gene</th>
                <th class="length-column">Length</th>
                <th class="non-atcg-column">Non-ATCG characters</th>
            </tr>
"""
            nonatcg_rows = ''.join([
                f"            <tr><td>{html.escape(r['gene'])}</td><td>{fmt_int(r['length'])}</td><td class='warnings'>{html.escape(r['non_atcg_chars'])}</td></tr>\n"
                for r in self.non_atcg_records
            ])
            nonatcg_section_html = nonatcg_table_header + nonatcg_rows + "        </table>\n"
        else:
            nonatcg_section_html = "<p style=\"color: #27ae60; font-weight: bold;\">✅ All sequences contain only canonical bases (ATCG)!</p>"

        # Codon Usage
        if self.participating_genes > 0:
            codon_usage_summary = (
                f"<div class='summary'><div class='summary-item'>Participating genes (ATCG-only): <strong>{fmt_int(self.participating_genes)}</strong></div>"
                f"<div class='summary-item'>Total codons (including STOP): <strong>{fmt_int(self.total_codon_count)}</strong></div>"
                f"<div class='summary-item'>Excluded genes (non-ATCG): <strong>{fmt_int(self.excluded_atcg_bad_genes)}</strong></div>"
                f"<div class='summary-item' title='Pre-filtered to ATCG-only, no skipped triplets generated'>Skipped triplets: <strong>0</strong></div></div>"
            )
            
            # Filters toolbar
            codon_usage_filters = """
        <div class="codon-usage-filters">
            <div class="filter-controls">
                <label for="aaFilter" class="filter-label">Amino acid:</label>
                <select id="aaFilter" class="filter-select" aria-label="Filter by amino acid">
                    <option value="">All</option>
                </select>
                
                <label for="searchFilter" class="filter-label">Search:</label>
                <input type="text" id="searchFilter" class="filter-input" placeholder="Search amino acid or codon..." aria-label="Search amino acid or codon">
                
                <button id="clearFilters" class="btn btn-secondary filter-btn" aria-label="Clear all filters">Clear filters</button>
                
                <span id="rowCount" class="row-count">Rows: 64 / 64</span>
            </div>
            
            <div class="export-controls">
                <label class="export-checkbox-label">
                    <input type="checkbox" id="exportFullDataset" class="export-checkbox">
                    Export full dataset
                </label>
                <button id="downloadCodonUsageXLSX" class="btn export-btn" aria-label="Export Codon Usage as XLSX with frozen header and AutoFilter">Export Codon Usage (XLSX)</button>
                <button id="downloadCodonUsageCSV" class="btn export-btn" aria-label="Export Codon Usage as CSV">Export Codon Usage (CSV)</button>
            </div>
            <div class="export-notes">
                <small>XLSX: header row is frozen and AutoFilter is enabled.</small><br>
                <small>CSV: plain text with header; filters and frozen header are not supported by CSV.</small>
            </div>
        </div>
"""
            
            codon_usage_table = """
        <div class="codon-usage-table-wrapper">
            <table id="codonUsageTable">
                <thead>
                    <tr>
                        <th class="codon-aa-column">Amino acid</th>
                        <th class="codon-column">Codon</th>
                        <th class="length-column">Count</th>
                        <th class="codon-freq-column">Frequency within AA (%)</th>
                    </tr>
                </thead>
                <tbody id="codonUsageTableBody">
                </tbody>
            </table>
        </div>
"""
            codon_usage_section_html = codon_usage_summary + codon_usage_filters + codon_usage_table
        else:
            codon_usage_section_html = "<p style=\"color: #e74c3c; font-weight: bold;\">⚠️ No compliant sequences available for codon usage</p>"

        # ---------- Main HTML (f-string; only inserts pre-built variables) ----------
        html_content = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
    <title>CDS Codon Usage Report and Codon Optimization</title>
    <!-- ExcelJS and FileSaver.js for XLSX export -->
    <script src="https://cdn.jsdelivr.net/npm/exceljs@4.4.0/dist/exceljs.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/file-saver@2.0.5/dist/FileSaver.min.js"></script>
    <style>{css_styles}</style>
</head>
<body>
    <div class=\"container\" id=\"pageTop\"> 
        <h1>CDS Codon Usage Report and Codon Optimization</h1>

        <nav class=\"report-nav\" id=\"reportNav\" aria-label=\"Report section navigation\">
            <a class=\"nav-link\" href=\"#summary\">Summary</a>
            <a class=\"nav-link\" href=\"#violations\">Violations</a>
            <a class=\"nav-link\" href=\"#invalid-records\">Invalid Records</a>
            <a class=\"nav-link\" href=\"#non-atcg\">Non-ATCG</a>
            <a class=\"nav-link\" href=\"#length-distribution\">Length Plot</a>
            <a class=\"nav-link\" href=\"#codon-usage\">Codon Usage</a>
            <a class=\"nav-link\" href=\"#upload-optimize\">Upload & Optimize</a>
            <a class=\"nav-link\" href=\"#methods-meta\">Methods & Meta</a>
        </nav>

        <!-- Summary -->
        <section id=\"summary\" class=\"report-section\">
        <h2>Summary</h2>
        <div class=\"summary\">
            <div class=\"summary-grid\">
                <div class=\"summary-card\">
                    <div class=\"summary-label\">Total records</div>
                    <div class=\"summary-value\">{fmt_int(self.total_records)}</div>
                </div>
                <div class=\"summary-card\">
                    <div class=\"summary-label\">Valid records</div>
                    <div class=\"summary-value valid\">{fmt_int(self.valid_records)}</div>
                </div>
                <div class=\"summary-card\">
                    <div class=\"summary-label\">Participating genes (ATCG-only)</div>
                    <div class=\"summary-value\">{fmt_int(self.participating_genes)}</div>
                </div>
                <div class=\"summary-card\">
                    <div class=\"summary-label\">Invalid records</div>
                    <div class=\"summary-value invalid\">{fmt_int(invalid_count)}</div>
                </div>
                <div class=\"summary-card\">
                    <div class=\"summary-label\">Invalid rate</div>
                    <div class=\"summary-value invalid\">{invalid_percentage:.2f}%</div>
                </div>
            </div>
            {"<div class='summary-warning'>⚠️ No sequence records parsed</div>" if self.total_records == 0 else ""}
        </div>
        </section>

        <!-- Violations Breakdown -->
        <section id=\"violations\" class=\"report-section\">
        <h2>Violations Breakdown</h2>
        <div class=\"violation-stats\">
            <div class=\"violation-item\">
                <div class=\"violation-number\">{fmt_int(self.violations['length'])}</div>
                <div class=\"violation-label\">Length violations</div>
                <div class=\"violation-rate\">Rate: {"—" if self.total_records == 0 else fmt_pct(self.violations['length']/self.total_records*100)}</div>
            </div>
            <div class=\"violation-item\">
                <div class=\"violation-number\">{fmt_int(self.violations['start_codon'])}</div>
                <div class=\"violation-label\">Start codon violations</div>
                <div class=\"violation-rate\">Rate: {"—" if self.total_records == 0 else fmt_pct(self.violations['start_codon']/self.total_records*100)}</div>
            </div>
            <div class=\"violation-item\">
                <div class=\"violation-number\">{fmt_int(self.violations['stop_codon'])}</div>
                <div class=\"violation-label\">Stop codon violations</div>
                <div class=\"violation-rate\">Rate: {"—" if self.total_records == 0 else fmt_pct(self.violations['stop_codon']/self.total_records*100)}</div>
            </div>
        </div>
        </section>

        <!-- Invalid Records -->
        <section id=\"invalid-records\" class=\"report-section\">
        <h2>Invalid Records</h2>
        {invalid_section_html}
        </section>

        <!-- Sequences with non-ATCG characters -->
        <section id=\"non-atcg\" class=\"report-section\">
        <h2>Sequences with non-ATCG characters</h2>
        {nonatcg_section_html}
        </section>

        <!-- Valid CDS Length Distribution -->
        {f"<section id='length-distribution' class='report-section'><h2>Valid CDS Length Distribution</h2><p class='section-note'>Distribution of CDS lengths among records that passed compliance checks.</p><div class='histogram-container'>{histogram_html}</div></section>" if histogram_html else ""}

        <!-- Codon Usage (including STOP) -->
        <section id=\"codon-usage\" class=\"report-section\">
        <h2>Codon Usage (including STOP)</h2>
        <p class=\"section-note\">Use filters to focus on specific amino acids or codons, then export the current view if needed.</p>
        {codon_usage_section_html}
        </section>

        <!-- Upload & Optimize -->
        <section id=\"upload-optimize\" class=\"report-section\">
        <h2>Upload & Optimize</h2>
        <div class="upload-card">
            <p>Only accepts a single-record CDS FASTA with .fa extension. The first line must start with '>' followed by the CDS name. Please upload a FASTA file ending with .fa.</p>
            <div class="controls">
                <input type="file" id="faUpload" accept=".fa" aria-label="Upload FASTA (.fa) file containing single CDS">
                <button id="btnOptimize" class="btn" disabled aria-label="Optimize based on the calculated codon usage">Optimize based on the calculated codon usage</button>
                <button id="btnClear" class="btn btn-secondary" aria-label="Clear upload and results">Clear</button>
            </div>
            <div id="statusToast" class="toast toast-info" role="status" aria-live="polite"></div>
            <div id="uploadResults" class="results"></div>
            <div id="optimizeResults" class="results"></div>
        </div>
        </section>

        <!-- Methods & Meta -->
        <section id=\"methods-meta\" class=\"report-section\">
        <h2>Methods & Meta</h2>
        <div class=\"meta\">
            <p><strong>Timestamp:</strong> {self.timestamp}</p>
            <p><strong>Author:</strong> Haoqiu Liu</p>
            <p><strong>Script version:</strong> {self.version}</p>
            <p><strong>Input file path:</strong> {os.path.abspath(self.input_file)}</p>
            <p><strong>Validation rules:</strong></p>
            <ul>
                <li>1. Sequence length must be a multiple of 3</li>
                <li>2. Start codon must be ATG</li>
                <li>3. Stop codon must be one of TGA, TAA, TAG (checked only when length is a multiple of 3)</li>
            </ul>
            <p><strong>Codon usage methodology:</strong></p>
            <ul>
                <li>Codon usage only includes ATCG-only compliant CDS; compliant CDS containing non-ATCG characters are completely excluded from codon usage; therefore Skipped triplets=0</li>
                <li>Statistical scope still includes start ATG and terminal STOP; synonymous codon grouping and STOP group remain unchanged</li>
                <li>Statistics count each 3bp triplet occurrence, including start ATG and stop TAA/TGA/TAG</li>
                <li>Synonymous codon aggregation follows standard genetic code mapping, with STOP treated as an independent amino acid group</li>
            </ul>
            <ul>
                <li>STOP is treated as an independent amino acid group for both statistics and optimization (replacement across STOP variants follows the top-frequency-within-group rule).</li>
                <li>Uploaded CDS is assumed to be on the coding strand (5'→3') and in the correct orientation; reverse-complement is not performed.</li>
            </ul>
            <p><strong>Note:</strong> Sequences containing non-ATCG characters are listed separately, but compliance determination still follows the above three rules.</p>
        </div>
        </section>

    </div>

    <button id="backToTopBtn" class="back-to-top" type="button" aria-label="Back to top">Top</button>

    <!-- Inject read-only data objects -->
    <script>
    window.CODON_TO_AA = Object.freeze({json.dumps(codon_to_aa_obj, ensure_ascii=False)});
    window.AA_CODON_FREQ = Object.freeze({json.dumps(aa_codon_freq_obj, ensure_ascii=False)});
    window.TOP_CODON_PER_AA = Object.freeze({json.dumps(top_codon_per_aa_obj, ensure_ascii=False)});
    window.AA_ONE_LETTER = Object.freeze({json.dumps(aa_one_letter_obj, ensure_ascii=False)});
    window.CODON_USAGE_ROWS = Object.freeze({json.dumps(codon_usage_rows_obj, ensure_ascii=False)});
    </script>

    <!-- Base scripts -->
    <script>{js_scripts}</script>
    
</body>
</html>
"""
        
        # Write HTML file to disk
        try:
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)
            print(f"HTML report generated: {os.path.abspath(output_file)}")
        except Exception as e:
            print(f"Error while generating report: {e}")
    
    def print_summary(self):
        """
        Print core statistical information to console
        """
        invalid_count = len(self.invalid_records)
        invalid_percentage = (invalid_count / self.total_records * 100) if self.total_records > 0 else 0
        
        print("\n" + "="*50)
        print("CDS Validation Summary")
        print("="*50)
        print(f"Total sequences: {self.total_records}")
        print(f"Valid sequences: {self.valid_records}")
        print(f"Invalid sequences: {invalid_count}")
        print(f"Invalid rate: {invalid_percentage:.2f}%")
        print(f"Sequences with non-ATCG characters: {len(self.non_atcg_records)}")
        print(f"Participating genes for codon usage (ATCG-only): {self.participating_genes}")
        print(f"Excluded genes (non-ATCG) for codon usage: {self.excluded_atcg_bad_genes}")
        print(f"Total codons (including STOP): {self.total_codon_count}")
        
        if self.total_records == 0:
            print("\nWarning: No sequence records parsed")
        else:
            print("\nViolation types:")
            print(f"  Length violations: {self.violations['length']} ({(self.violations['length']/self.total_records*100):.2f}%)")
            print(f"  Start codon violations: {self.violations['start_codon']} ({(self.violations['start_codon']/self.total_records*100):.2f}%)")
            print(f"  Stop codon violations: {self.violations['stop_codon']} ({(self.violations['stop_codon']/self.total_records*100):.2f}%)")
        print("="*50)


def main():
    """
    Main function
    """
    print("CDS Codon Usage Checker v1.0")
    print("Starting...")
    
    start_time = time.time()
    
    # Create checker instance
    checker = CDSChecker()
    
    # Check if input file exists
    if not os.path.exists(checker.input_file):
        print(f"❌ Input file not found: {checker.input_file}")
        return
    
    try:
        # Process sequences
        checker.process_sequences()
        
        # Generate HTML report even if no records parsed
        checker.generate_html_report()
        
        # Print summary
        checker.print_summary()
        
        end_time = time.time()
        print(f"\n✅ Done. Elapsed time: {end_time - start_time:.2f} s")
        
    except KeyboardInterrupt:
        print("\nInterrupted by user")
    except Exception as e:
        print(f"Error during execution: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
