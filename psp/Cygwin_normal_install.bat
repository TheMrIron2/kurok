@echo off

D:
chdir D:\cygwin\bin

bash --login -i -c "cd \"Kurok\psp\"; make -f MakefileNormal install"
