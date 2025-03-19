# python_cli_template
 My template for python cli tools


Explanations in the scripts.

# Install

Using `pip` to install any python "package" as a command line tool:

This requires the `pyproject.toml` file included in this repo. Please see notes therein. Setuptools may not be the best package manager nowadays. Please feel free to explore others.

From the terminal:

`cd python_cli_template`

`pip install .`

* This should install `pct` as a runnable command from the terminal

Run these on the command line to get started.

"DEV" - this is how you'll run things when you're making updates and changes (specify path to script)

"INSTALLED" - this is how you'll run things after `pip` installing a stable version (from any directory)

### Basic

DEV:

`python src/pct/pct.py`

INSTALLED:

`pct`

### Help Menu

DEV:

`python src/pct/pct.py -h`

INSTALLED:

`pct -h`

### `shuth` subcommand (help menu)

DEV:

`python src/pct/pct.py shuth -h`

INSTALLED:

`pct shuth -h`


### An arbitrary example of actually running something:

DEV:

`python src/pct/pct.py sqawk tigers -c title`

INSTALLED:

`pct sqawk tigers -c title`


returns:

`quietly: Tigers!!!!`