#Reading CSV file using Pandas

import pandas as pd
import os 
from Bio import SeqIO
import matplotlib.pyplot as plt

base_dir = "/Users/admin/Desktop/Repositories/python_for_biologists/datasets/" 
file_name = 'all_weather.csv'
file_path = os.path.join(base_dir, file_name)
print(f'File path: {file_path}')

if os.path.exists(file_path):
    print('File found. Reading...')
    df = pd.read_csv(file_path)
    print(df.head())
else:
    print("Error: File not found at", file_path)
    

#Parse/Analyze FASTA file. Count total sequences, print length of each seq.

cwd = os.getcwd()
print('cwd:', cwd)
new_dir = '/Users/admin/Desktop/Repositories/python_for_biologists/BRCA2_datasets/ncbi_dataset/data'
os.chdir(new_dir)
print('pwd:', new_dir)

fasta_file = "gene.fasta"

for record in SeqIO.parse(fasta_file, 'fasta'):
    sequence = record.seq
    
    print(f'ID: {record.id}')
    print(f'Description: {record.description}')
    print(f'Sequence: {record.seq}')
    print(f'Length: {len(record.seq)}\n')

#Write a script to calculate GC content of a sequence

def gc_content (sequence):
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100

gc_percentage = print(f'The GC % content: {gc_content(sequence)}'+'%')

#Visualize sample data using matplotlib

sequence_lengths = []
sequence_ids = []

for record in SeqIO.parse(fasta_file, "fasta"):
    sequence_ids.append(record.id)        # Sequence ID
    sequence_lengths.append(len(record.seq))  # Length of the sequence


plt.figure(figsize=(10, 6))
plt.plot(sequence_ids, sequence_lengths, marker='o', linestyle='-', color='b', label='Sequence Length')

plt.xlabel('Sequence ID')
plt.ylabel('Sequence Length')
plt.title('Sequence Lengths of Sequences in FASTA File')
plt.xticks(rotation=45, ha='right')  # Rotate the x-axis labels for better visibility

plt.tight_layout()  # Adjust layout for better spacing
plt.legend()
plt.show()
