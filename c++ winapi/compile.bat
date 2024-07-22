if exist *.exe del *.exe
if exist *.dll del *.dll
cl rgz_vatlin.cpp kernel32.lib user32.lib gdi32.lib advapi32.lib
cl /LD info.cpp kernel32.lib user32.lib gdi32.lib advapi32.lib
del *.obj *.lib *.exp