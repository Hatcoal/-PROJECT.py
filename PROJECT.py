import tkinter as tk
from tkinter import filedialog, simpledialog

def translate_DNA(DNA_seq):
    codon_table = {
        'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
        'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
        'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stop', 'TAG': 'Stop',
        'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stop', 'TGG': 'Trp',
        'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
        'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
        'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
        'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    triplets_frame1 = []
    triplets_frame2 = []
    triplets_frame3 = []

    protein_seq_frame1 = []
    protein_seq_frame2 = []
    protein_seq_frame3 = []

    for i in range(len(DNA_seq)):
        triplet = DNA_seq[i:i + 3]
        if len(triplet) == 3:
            amino_acid = codon_table.get(triplet, "X")

            if i % 3 == 0:
                protein_seq_frame1.append(amino_acid)
                triplets_frame1.append(triplet)
            elif (i - 1) % 3 == 0:
                protein_seq_frame2.append(amino_acid)
                triplets_frame2.append(triplet)
            elif (i - 2) % 3 == 0:
                protein_seq_frame3.append(amino_acid)
                triplets_frame3.append(triplet)

    spaced_seq_frame1 = ' '.join(protein_seq_frame1)
    spaced_seq_frame2 = ' '.join(protein_seq_frame2)
    spaced_seq_frame3 = ' '.join(protein_seq_frame3)
    spaced_triplets_frame1 = ' '.join(triplets_frame1)
    spaced_triplets_frame2 = ' '.join(triplets_frame2)
    spaced_triplets_frame3 = ' '.join(triplets_frame3)

    return (spaced_seq_frame1, spaced_triplets_frame1,
            spaced_seq_frame2, spaced_triplets_frame2,
            spaced_seq_frame3, spaced_triplets_frame3)

def reverse_complement(DNA_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    complement_seq = ''.join(complement[base] for base in reversed(DNA_seq))

    return complement_seq

def translate_reverse_complement(DNA_seq):
    complement_seq = reverse_complement(DNA_seq)

    (Protein_seq_frame1, Triplets_seq1,
     Protein_seq_frame2, Triplets_seq2,
     Protein_seq_frame3, Triplets_seq3) = translate_DNA(complement_seq)

    result_text.config(state="normal")
    result_text.delete("1.0", "end")

    result_text.insert("end", "\nReverse Complement:\n")
    result_text.insert("end", complement_seq + "\n")

    result_text.insert("end", "\nFrame 1:    " + Triplets_seq1 + "\nProtein Seq: " + Protein_seq_frame1 + "\n")
    result_text.insert("end", "\nFrame 2:    " + Triplets_seq2 + "\nProtein Seq: " + Protein_seq_frame2 + "\n")
    result_text.insert("end", "\nFrame 3:    " + Triplets_seq3 + "\nProtein Seq: " + Protein_seq_frame3 + "\n")

    result_text.config(state="disabled")

def find_start_stop_codon(dna_seq, start_codon='ATG', stop_codons=['TAA', 'TAG', 'TGA']):
    codon_table = {
        'TTT': 'Phe', 'TTC': 'Phe', 'TTA': 'Leu', 'TTG': 'Leu',
        'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
        'TAT': 'Tyr', 'TAC': 'Tyr', 'TAA': 'Stop', 'TAG': 'Stop',
        'TGT': 'Cys', 'TGC': 'Cys', 'TGA': 'Stop', 'TGG': 'Trp',
        'CTT': 'Leu', 'CTC': 'Leu', 'CTA': 'Leu', 'CTG': 'Leu',
        'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
        'CAT': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
        'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
        'ATT': 'Ile', 'ATC': 'Ile', 'ATA': 'Ile', 'ATG': 'Met',
        'ACT': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAT': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
        'AGT': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
        'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val', 'GTG': 'Val',
        'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
        'GAT': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
        'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
    }

    translations = []
    in_frame = False

    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i + 3]

        if codon == start_codon:
            in_frame = True
            translations.append('Met')
        elif codon in stop_codons and in_frame:
            break
        elif in_frame:
            amino_acid = codon_table.get(codon, "X")
            translations.append(amino_acid)

    return ' '.join(translations)

def on_open_file():
    file_path = filedialog.askopenfilename(title="Select a DNA Sequence File", filetypes=(("Fasta Files", "*.fasta"), ("All Files", "*.*")))
    if file_path:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            dna_sequence = "".join(line.strip() for line in lines[1:]) if len(lines) > 1 else lines[0].strip()
            text_box.delete("1.0", "end")
            text_box.insert("1.0", dna_sequence)


def on_clear():
    text_box.delete("1.0", "end")

def on_reset():
    text_box.delete("1.0", "end")
    result_text.config(state="normal")
    result_text.delete("1.0", "end")
    result_text.config(state="disabled")

def save_results():
    result = result_text.get("1.0", "end-1c")
    if result:
        file_path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text Files", "*.txt")])
        if file_path:
            with open(file_path, 'w') as file:
                file.write(result)

def on_translate():
    sequence_type = simpledialog.askstring("Input", "Are you providing DNA or RNA sequence? (Type 'DNA' or 'RNA')")

    if sequence_type and sequence_type.upper() in ['DNA', 'RNA']:
        DNA_Seq = text_box.get("1.0", "end-1c").upper()

        result_text.config(state="normal")
        result_text.delete("1.0", "end")

        if all(base in 'ATCG' for base in DNA_Seq):
            result_text.insert("end", f"{sequence_type} Sequence is valid.\n")

            start_stop_translation = find_start_stop_codon(DNA_Seq)
            result_text.insert("end", f"{sequence_type} Seq: {DNA_Seq}\n")

            if start_stop_translation:
                result_text.insert("end", f"\nTranslation Between Start and Stop Codons:\n{start_stop_translation}\n")
            else:
                result_text.insert("end", "No start and stop codon pairs found.\n")

            if sequence_type.upper() == 'DNA':
                translate_reverse_complement(DNA_Seq)
            else:
                RNA_Seq = DNA_Seq.replace('T', 'U')
                translate_reverse_complement(RNA_Seq)
        else:
            result_text.insert("end", f"Error: Invalid {sequence_type} sequence.")

app = tk.Tk()
app.title("FA20-TSU-777| JOSHUA BENSON")
# app.geometry("800x600")
app.configure(bg="skyblue")
app.minsize(height=700,width=700)
app.maxsize(height=700,width=700)
frame_input = tk.Frame(app, padx=10, pady=10)
frame_input.grid(row=0, column=0, sticky="nsew")

frame_buttons = tk.Frame(app)
frame_buttons.grid(row=1, column=0, pady=10, sticky="nsew")

frame_results = tk.Frame(app, padx=10, pady=10)
frame_results.grid(row=2, column=0, sticky="nsew")

text_box = tk.Text(frame_input, height=6, width=60, font=("Arial", 12))
text_box.grid(row=0, column=0, padx=10, pady=10)
text_box.insert("1.0", "Enter your DNA/RNA sequence here...")

open_file_button = tk.Button(frame_buttons, text="Open File", command=on_open_file, bg="lightgreen", font=("Arial", 12), padx=10, pady=5)
open_file_button.grid(row=0, column=0, padx=10, pady=10)

translate_button = tk.Button(frame_buttons, text="Translate", command=on_translate, bg="lightblue", font=("Arial", 12), padx=10, pady=5)
translate_button.grid(row=0, column=1, padx=10, pady=10)

clear_button = tk.Button(frame_buttons, text="Clear", command=on_clear, bg="lightcoral", font=("Arial", 12), padx=10, pady=5)
clear_button.grid(row=0, column=2, padx=10, pady=10)

reset_button = tk.Button(frame_buttons, text="Reset", command=on_reset, bg="lightyellow", font=("Arial", 12), padx=10, pady=5)
reset_button.grid(row=0, column=3, padx=10, pady=10)

result_text = tk.Text(frame_results, height=12, width=80, state="disabled", font=("Arial", 12), bg="lightgray")
result_text.grid(row=0, column=0, pady=(0, 10))

save_button = tk.Button(frame_results, text="Save Results", command=save_results, bg="lightgreen", font=("Arial", 12), padx=10, pady=5)
save_button.grid(row=1, column=0, pady=(0, 10))

instructions_label = tk.Label(app, text="Enter your DNA/RNA sequence or open a file.", font=("Arial", 12))
instructions_label.grid(row=3, column=0, pady=(10, 0))

app.update_idletasks()
app.geometry(f"{app.winfo_width()}x{app.winfo_height()}+{app.winfo_screenwidth() // 2 - app.winfo_width() // 2}+{app.winfo_screenheight() // 2 - app.winfo_height() // 2}")

app.mainloop()
