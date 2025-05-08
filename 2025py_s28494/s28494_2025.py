#!/usr/bin/env python3
"""
FASTA DNA Sequence Generator

Ten skrypt generuje losowe sekwencje DNA (A, C, G, T) w formacie FASTA.

Cel programu:
- Demonstracja generacji sekwencji DNA z interaktywną konfiguracją,
  w tym niestandardowym tagiem identyfikacyjnym, walidacją ID,
  zapisem do katalogu i podstawowymi statystykami.
Kontekst zastosowania:
- Bioinformatyka edukacyjna, szybkie prototypowanie danych DNA
  lub testowanie pipeline'ów analiz sekwencji.
"""

# Import bibliotek standardowych
import random      # do losowania nukleotydów i pozycji
import argparse    # do obsługi argumentów wiersza poleceń
import os          # do obsługi operacji na ścieżkach i katalogach

# === Ulepszenie 1: niestandardowa nazwa tagu ===
# ORIGINAL:
# TAG_NAME = "AAA"
# MODIFIED (pozwala użytkownikowi ustawić własny tag, domyślnie 'AAA'):
TAG_NAME = input('Podaj nazwę tagu (domyślnie AAA): ').strip() or 'AAA'
# Dlaczego: większa elastyczność i możliwość brandingu generowanych plików

# Stałe formatujące
INSERT_WIDTH = len(TAG_NAME)   # długość tagu, używana przy wstawianiu
LINE_WIDTH = 60               # szerokość linii w pliku FASTA


def generate_dna_sequence(length: int) -> str:
    """
    Generuje losową sekwencję DNA o podanej długości.
    :param length: liczba nukleotydów do wygenerowania
    :return: ciąg znaków złożony z A, C, G i T
    """
    return ''.join(random.choices('ACGT', k=length))


def compute_stats(sequence: str) -> dict:
    """
    Oblicza statystyki nukleotydów dla podanej sekwencji (bez tagu).
    :param sequence: oryginalna sekwencja DNA
    :return: słownik z liczbami, procentami i stosunkiem C+G / A+T
    """
    length = len(sequence)
    # Zliczanie każdego nukleotydu
    counts = {nuc: sequence.count(nuc) for nuc in 'ACGT'}
    # Procentowa zawartość każdego nukleotydu
    percents = {nuc: (counts[nuc] / length) * 100 for nuc in counts}
    cg = counts['C'] + counts['G']
    at = counts['A'] + counts['T']
    ratio = (cg / at) if at > 0 else None
    return {'counts': counts, 'percents': percents, 'cg_at_ratio': ratio}


def format_fasta(seq_id: str, sequence: str, description: str = '') -> str:
    """
    Formatuje sekwencję wraz z tagiem do formatu FASTA.
    :param seq_id: identyfikator sekwencji
    :param sequence: sekwencja DNA z wstawionym tagiem
    :param description: tekst opisu umieszczany w nagłówku FASTA
    :return: gotowy tekst w formacie FASTA
    """
    header = f">{seq_id} {description}".strip()
    # Zawijanie sekwencji co LINE_WIDTH znaków
    wrapped = '\n'.join(
        sequence[i:i+LINE_WIDTH] for i in range(0, len(sequence), LINE_WIDTH)
    )
    return f"{header}\n{wrapped}\n"


def main():
    # Pobranie długości sekwencji przez input
    try:
        length = int(input('Podaj długość sekwencji (liczba całkowita > 0): ').strip())
        if length <= 0:
            raise ValueError
    except ValueError:
        print('Nieprawidłowa wartość długości. Podaj dodatnią liczbę całkowitą.')
        return

    # === Ulepszenie 2: obsługa katalogu wyjściowego ===
    parser = argparse.ArgumentParser(
        description='Generuj losowe sekwencje DNA i zapisuj je do plików FASTA.'
    )
    # ORIGINAL:
    # parser.add_argument('-n', '--number', type=int, default=1,
    #                     help='Liczba sekwencji do wygenerowania')
    # MODIFIED (dodano opcję katalogu wyjściowego oraz przenumerowano argumenty):
    parser.add_argument('-n', '--number', type=int, default=1,
                        help='Liczba sekwencji do wygenerowania')
    parser.add_argument('-d', '--outdir', default='.',
                        help='Katalog wyjściowy dla plików FASTA')
    args = parser.parse_args()

    # Tworzenie katalogu, jeśli nie istnieje
    os.makedirs(args.outdir, exist_ok=True)

    for i in range(1, args.number + 1):
        # Pobranie i walidacja ID sekwencji
        # === Ulepszenie 3: walidacja ID ===
        # ORIGINAL:
        # seq_id = input(f'Podaj ID dla sekwencji {i}: ').strip()
        # MODIFIED (ID musi zawierać tylko litery, cyfry lub podkreślenia):
        raw_id = input(f'Podaj ID dla sekwencji {i}: ').strip()
        if not raw_id or not raw_id.replace('_', '').isalnum():
            print('Nieprawidłowe ID. Użyj tylko liter, cyfr i podkreśleń.')
            return
        seq_id = raw_id

        # Pobranie opisu
        description = input(f'Podaj opis dla sekwencji {i} (opcjonalnie): ').strip()

        # Generacja oryginalnej sekwencji bez tagu
        seq = generate_dna_sequence(length)
        stats = compute_stats(seq)  # statystyki na oryginale

        # Wstawienie tagu w losowym miejscu
        insert_pos = random.randint(0, len(seq))
        seq_with_tag = seq[:insert_pos] + TAG_NAME + seq[insert_pos:]

        # Formatowanie i zapis do pliku
        fasta_content = format_fasta(seq_id, seq_with_tag, description)
        # Ścieżka pliku: katalog + nazwa ID.fasta
        filepath = os.path.join(args.outdir, f"{seq_id}.fasta")
        try:
            with open(filepath, 'w') as f:
                f.write(fasta_content)
            print(f"Zapisano sekwencję do pliku: {filepath}")
        except IOError as e:
            print(f"Błąd zapisu pliku {filepath}: {e}")
            continue

        # Wyświetlenie statystyk (bez tagu)
        print(f"\nStatystyki dla sekwencji {seq_id} (bez tagu '{TAG_NAME}'):")
        for nuc, pct in stats['percents'].items():
            print(f"  {nuc}: {pct:.2f}% ({stats['counts'][nuc]} nukleotydów)")
        ratio = stats['cg_at_ratio']
        if ratio is not None:
            print(f"  Stosunek (C+G)/(A+T): {ratio:.2f}\n")
        else:
            print("  Stosunek (C+G)/(A+T): nieokreślony (brak A lub T)\n")


if __name__ == '__main__':
    # Uruchomienie programu
    main()