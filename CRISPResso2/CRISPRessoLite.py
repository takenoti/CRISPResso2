from CRISPResso2 import CRISPResso2Align, alignment_tags
from typing_extensions import Annotated
import functools
import numpy
import pathlib
import sys
import typer


def process_reads(
    amplicon_sequence: str,
    cut_site_offset: Annotated[
        int,
        typer.Argument(
            help="zero-based position of the nucleotide left to the cutsite"
        ),
    ] = None,
    input_file: Annotated[typer.FileText, typer.Option()] = sys.stdin,
    output_file: Annotated[typer.FileTextWrite, typer.Option()] = sys.stdout,
    gap_incentive_score: Annotated[
        str,
        typer.Option(
            help="a scalar gap_incentive_score, or a vector of length len(amplicon)+1 with per-position gap incentive scores."
            "position 0 refers to gaps prior to the amplicon start (zero-based indexing)."
        ),
    ] = "1",
    gap_open_score: int = -20,
    gap_extend_score: int = -2,
    alignment_matrix_path: pathlib.Path = pathlib.Path(__file__).parent / "EDNAFULL",
    short_cs_tag: bool = False,
):
    alignment_matrix = CRISPResso2Align.read_matrix(alignment_matrix_path)

    # Create gap incentive around the cutsite
    if gap_incentive_score.isnumeric():
        gap_incentive = numpy.array(
            [0] * (cut_site_offset)
            + [int(gap_incentive_score)]
            + [0] * (len(amplicon_sequence) - cut_site_offset)
        )
    else:
        if cut_site_offset is not None:
            raise ValueError(
                "Cannot specify both 'cut_site_offset' and a vector of gap incentive scores!"
            )
        gap_incentive = numpy.array([int(i) for i in gap_incentive_score.split(",")])
        if len(gap_incentive) != len(amplicon_sequence) + 1:
            raise ValueError(
                f"Length of gap incentive score {len(gap_incentive)} does not match amplicon length {len(amplicon_sequence)}! "
                "gap_incentive_score must be longer than amplicon by exactly 1bp. "
                f"gap_incentive_score={gap_incentive_score}, amplicon_sequence='{amplicon_sequence}'"
            )

    @functools.lru_cache()
    def process_read_sequence(
        amplicon_sequence: str,
        cut_site_offset: int,  # zero-based
        read_sequence: str,
    ):
        aligned_read_sequence, aligned_amplicon, alignment_score = (
            CRISPResso2Align.global_align(
                read_sequence,
                amplicon_sequence,
                matrix=alignment_matrix,
                gap_incentive=gap_incentive,
                gap_open=gap_open_score,
                gap_extend=gap_extend_score,
            )
        )

        return (
            alignment_tags.aligned_crispresso_to_cs_tag(
                aligned_amplicon,
                aligned_read_sequence,
                short_cs_tag=short_cs_tag,
            ),
            alignment_score,
        )

    # read in fastq format
    while True:
        header = input_file.readline()
        if not header:
            break

        if header[0] == ">":
            # FASTA
            read_id = header[1:-1]  # Remove `>` and `\n`
            read_sequence = input_file.readline().rstrip()
        else:
            # FASTQ
            read_id = header[1:-1].split(" ", 1)[0]
            read_sequence = input_file.readline().rstrip()
            _plus = input_file.readline()
            _quality = input_file.readline()

        cs_tag, alignment_score = process_read_sequence(
            amplicon_sequence=amplicon_sequence,
            cut_site_offset=cut_site_offset,
            read_sequence=read_sequence,
        )

        output_file.write(
            ",".join(
                (
                    read_id,
                    cs_tag,
                    str(alignment_score),
                    "\n",
                )
            )
        )


def main():
    typer.run(process_reads)


if __name__ == "__main__":
    main()
