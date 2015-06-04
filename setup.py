import os
from distutils.core import setup
import py2exe

print "Removing Trash"
os.system("rmdir /s /q build")

# delete the old build drive
os.system("rmdir /s /q dist")


setup(
    options={
        "py2exe":{
            "dll_excludes": ["MSVCP90.dll", "HID.DLL", "w9xpopen.exe",
                'libifcoremd.dll', ],#'msvcr71.dll'
            #"ascii": 1, # to make a smaller executable, don't include the encodings
            "compressed": 0, # compress the library archive
            #"skip_archive": 1, 
        }
    },
    console=[{'script': 'baikal2gse.py'}],
    data_files=[('.', ['libgse2.pyd'])],
    includes=['UserDict', 'UserList', 'UserString', 'warnings', '__future__',
        'future', 'future.builtins'],
)
