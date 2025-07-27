def parse_fna(filepath):
    '''
    Description: 
    Parses a given .fna file and returns its indexable contents

    IN:
    para1: path to .fna file; [string]

    OUT:
    return1: sequence data [(sequence_id, sequence_size, sequence)]; [list[tuple(3)]]
            where sequence_id is in format "<SequenceName SequenceNumber"
    '''
    sequences = []
    seq_id = None
    seq_lines = []

    with open(filepath, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # if weve already seen the sequence id or are mid sequence, parse the rest
                if seq_id is not None:
                    full_seq = ''.join(seq_lines)
                    sequences.append((seq_id, len(full_seq), full_seq))
                # encountering a new sequence id
                seq_id = line.strip()  # strip whitespace
                seq_lines = []
            else:
                seq_lines.append(line)

        # the last sequence
        if seq_id is not None:
            full_seq = ''.join(seq_lines)
            sequences.append((seq_id, len(full_seq), full_seq))

    return sequences


def write_pairs(output_file, seq1_data, seq2_data, score):
    '''
    Description:
    Writes two aligned sequences to a .fna file in standard FASTA format.

    IN:
    seq1_data: tuple containing (sequence_id, aligned_sequence) [tuple(str, str)]
            sequence_id should be in the format '>SequenceName SequenceNumber'
    seq2_data: tuple containing (sequence_id, aligned_sequence) [tuple(str, str)]
    score: alignment score [int]

    OUT:
    Writes a .fna file named 'aligned_output.fna' in the current directory.
    '''

    sequence_1_id, sequence_1_aligned = seq1_data
    sequence_2_id, sequence_2_aligned = seq2_data

    with open(output_file, "w") as f:
        f.write(f"{sequence_1_id+'; score='+str(score)}\n")
        f.write(f"{sequence_1_aligned}\n")
        f.write(f"{sequence_2_id+'; score='+str(score)}\n")
        f.write(f"{sequence_2_aligned}\n")


def main():
    fna_path = './src/input/7.fna'
    output_path = './data/parsefnadebug.txt'

    sequence_data = parse_fna(fna_path)
    num_sequences = len(sequence_data)

    print(f'Success')
    with open(output_path, 'w', encoding='utf-8') as file:
        file.write('Sequence Data : '+str(num_sequences)+' Sequences\n')
        file.write(str(sequence_data) + '\n')


if __name__ == "__main__":
    main()
