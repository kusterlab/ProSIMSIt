import tomli
import logging
import argparse
from pathlib import Path

from . import __version__, __copyright__

logger = logging.getLogger(__name__)


class ArgumentParserWithLogger(argparse.ArgumentParser):
    def error(self, message):
        logger.error(f"Error parsing input arguments: {message}")
        super().error(message)


def parse_args(argv):
    desc = f'ProSIMSIt version {__version__}\n{__copyright__}'
    apars = ArgumentParserWithLogger(
        description=desc, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument(
        "-c",
        "--config_path",
        default=None,
        metavar="DIR",
        required=True,
        help=(
            "Path to config file in json format."
        ),
    )
    args = apars.parse_args()
    return args


def read_config(argv):
    config_path = parse_args(argv).config_path

    logger.info(f"Reading configuration from {config_path}")
    if isinstance(config_path, str):
        config_path = Path(config_path)
    with open(config_path, mode='rb') as f:
        data = tomli.load(f)
    return data


if __name__ == '__main__':
    raise NotImplementedError('Do not run this script.')
