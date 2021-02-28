#!/usr/bin/env python

from .helper import get_date

def main(argv=None):
    """
    Function for print the date         
    """

    today = get_date()
    print_statement = f"Hello world. Today is {today}"
    print(print_statement)

    
if __name__ == "__main__":
    main(sys.argv[1:])
