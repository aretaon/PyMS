# -*- coding: utf-8 -*-
"""
Methods related to string manipulation

@author: aretaon
"""


def alphanum_string(s):
    import re
    # new compiler that finds non-alphanumeric characters
    rex = re.compile(r'\W')
    # actually replace the strings
    result = rex.sub('_', s)

    # remove double occurences of _ in the string
    # initialize
    old_char = ''
    new_result = ''
    for char in result:
        if old_char != char:
            new_result += char
        old_char = char

    return new_result
