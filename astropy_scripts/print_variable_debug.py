#Use in python (script) as:
# from __future__ import (absolute_import, division, print_function,
#                         unicode_literals) # python 2 as python 3
# from print_variable_debug import print_variable_debug as var_debug
# print("some_var")
# var_debug(some_var)

from __future__ import (absolute_import, division, print_function,
                        unicode_literals) # python 2 as python 3

def print_variable_debug(var, DEBUG=False):
    """Print debug info for a specific variable.
    """
    # TODO: var name printing is not workin!!!
    # refs
    # unknown var name: http://stackoverflow.com/questions/2273211/how-to-define-generic-variables-in-python-syntax-question
    # genergic var type: http://stackoverflow.com/questions/1515412/declaring-unknown-type-variable-in-python
    #
##    ##var = type(var)()
##    # print var name:
##    for k, v in list(locals().iteritems()):
##        if id(v) == id(var):
##            var_as_str = k
##
##    print(var_as_str)
##    #print(type(var_as_str)) # type of var name
##
    print(var)
    print(type(var))
    #print(str(var)) # equivalent to print(var)
    print(repr(var)) #nice dict, as if calling "var" in python interactive seesion (ipython)

    # detect if the var has a dict associated
    no_dict = False
    #no_dict_var_types = ['float', 'dict']
    #if type(var) in no_dict_var_types: ##not working!!!!
    if type(var) is float or type(var) is dict:
        no_dict = True

    if not no_dict:
        print(var.__dict__) #this is ugly! I should use it onlyin DEBUG, but the nice one should work 1st!!! #NOT ALL TYPES HAVE A DICT!!!! (i.e. float)
    if DEBUG:
        print("id =", id(var)) #memory address
