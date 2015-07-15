#!bin/bash/

file="LEER"

a=".tex"
b=".dvi"
c=".ps"

latex $file$a # crea el archivo file.dvi

dvips -Ppdf -G0 $file$b # crea el archio file.ps

ps2pdf $file$c # crea el archivo file.pdf
