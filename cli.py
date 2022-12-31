import click
import sys
from irms.data_analysis.data_analysis import label_aa_peaks


@click.group(
    name="program",
    subcommand_metavar="COMMAND <args>",
    short_help="Does something with something",
    context_settings=dict(max_content_width=85, help_option_names=["-h", "--help"]),
)
def program_cli():
    print("whoooo we did it")



@click.command(
    name="command",
    short_help="short help",
    help="long help",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
@click.option(
    "-a",
    "--add-labels",
    "some_function",
    help="When present it will add do something",
    required=False,
    default=None,
    type=click.Path(exists=False, dir_okay=False, writable=True),
)
@click.option(
    "-i", "--ignore", help="Ignore something", is_flag=True
)
@click.option(
    "-v", "--verbose", help="When present will set logging level to debug", is_flag=True
)
@click.option(
    "-f", "--fail", help="When present will set tell program to fail", is_flag=True, default=False
)
def command(input_file, output_file, fail):
    print("running command")
    if fail:
        sys.exit(1)
    else:
        sys.exit(0)

@click.command(
    name="analyze_ms_data",
    short_help="short help",
    help="long help",
)
@click.argument("input_file", nargs=1, type=click.Path(exists=True, dir_okay=False))
@click.argument("output_file", nargs=1, type=click.Path(exists=False, dir_okay=False))
def analyze_ms_data(input_file, output_file):
    label_aa_peaks(input_file, write=True, output_file=output_file)


program_cli.add_command(command)
program_cli.add_command(analyze_ms_data)


if __name__ == "__main__":
    program_cli()