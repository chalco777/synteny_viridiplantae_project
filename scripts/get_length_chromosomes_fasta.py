#!/usr/bin/env python3

"""
Script Function: Calculate the lengths of sequences in a FASTA file and output them in a TSV file, sorted by descending length.
"""

import os
import sys

def calcular_longitudes(fasta_path, output_path):
    """
    Reads a FASTA file, calculates the lengths of all sequences, and outputs them to a TSV file.

    Parameters:
        fasta_path (str): Path to the input FASTA file.
        output_path (str): Path to the output TSV file.

    Functionality:
        - Parses the FASTA file to extract sequence lengths.
        - Sorts sequences by length in descending order.
        - Writes the sequence names and lengths to a TSV file.
    """
    if not os.path.isfile(fasta_path):
        print(f"Error: The FASTA file '{fasta_path}' does not exist.")
        sys.exit(1)

    secuencias = {}
    with open(fasta_path, 'r') as f:
        nombre_secuencia = ''
        secuencia = ''

        for linea in f:
            linea = linea.strip()
            if linea.startswith('>'):  # New sequence header
                if nombre_secuencia:
                    secuencias[nombre_secuencia] = len(secuencia)
                # Extract the sequence name (identifier only, no extra info)
                nombre_secuencia = linea[1:].split()[0]
                secuencia = '' 
            else:
                secuencia += linea  

        # Store the last sequence
        if nombre_secuencia:
            secuencias[nombre_secuencia] = len(secuencia)

    # Sort sequences by length in descending order
    secuencias_ordenadas = sorted(secuencias.items(), key=lambda x: x[1], reverse=True)

    # Write the output TSV file
    with open(output_path, 'w') as out_f:
        out_f.write("Contig\tLongitud\n")
        for nombre, longitud in secuencias_ordenadas:
            out_f.write(f"{nombre}\t{longitud}\n")

    print(f"TSV file successfully generated: {output_path}")

if __name__ == "__main__":
    """
    Main script execution:
        - Calculates sequence lengths from the given FASTA file.
        - Outputs a sorted TSV file with sequence lengths.
    """
    if len(sys.argv) != 3:
        print("Usage: python calcular_longitudes.py /absolute/path/input.fasta /absolute/path/output.tsv")
        sys.exit(1)

    ruta_fasta = sys.argv[1]
    ruta_salida = sys.argv[2]

    calcular_longitudes(ruta_fasta, ruta_salida)
