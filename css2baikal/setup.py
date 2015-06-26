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
            "dll_excludes": ["MSVCP90.dll", "MSVCM90.DLL", "MSVCR90.DLL",
                "HID.DLL", "w9xpopen.exe", 'libifcoremd.dll'],
            "compressed": 0, # compress the library archive
            "excludes":['IPython', 'Image', 'PIL', 'Tkconstants', 'Tkinter',
                '_hashlib', '_imaging', '_ssl', '_ssl', 'bz2', 'calendar',
                'compiler', 'cookielib', 'cookielib', 'doctest',
                'email', 'nose', 'optparse', 'pdb', 'pydoc', 'pywin',
                'readline', 'sqlite3', 'tcl', 'tornado', 'zipfile', 'zmq'],
        }
    },
    console=[{'script': 'css2baikal.py'}],
    data_files=[('.', ['multiarray.pyd'])],
)
