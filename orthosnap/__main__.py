"""orthosnap.__main__: executed when orthosnap is called as script"""

import sys

from .orthosnap import main

main(sys.argv[1:])
