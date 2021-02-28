#!/usr/bin/env python

from .helper import get_date

def main(argv=None):
    """
    Function for print the date
    """

    date_and_time = get_date()
    print_statement = f"Hello world. The date and time is {date_and_time}."
    print(print_statement)

if __name__ == "__main__":
    main(sys.argv[1:])
