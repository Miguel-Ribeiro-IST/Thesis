REM -----------------------------------------------------------------
REM Compilation on TeXworks
REM
REM script to be used only when the glossary or nomenclature are modified
REM otherwise, command "texify" inside TeXworks does the job better
REM
REM the defaul is the optimized full compilation option
REM uncommented other command options as desired
REM
REM Alexandre Gomes 4legom@gmail.com
REM 28 September 2015
REM -----------------------------------------------------------------

REM -----------------------------------------------------------------
REM Optimized full compilation

@echo off

echo. &
texify --pdf Thesis.tex
echo. &
makeindex Thesis.glo -s Thesis.ist -t Thesis.glg -o Thesis.gls
echo. &
makeindex Thesis.nlo -s nomencl.ist -o Thesis.nls
echo. &
texify --pdf Thesis.tex


REM -----------------------------------------------------------------
REM Glossary and nomenclature only
REM 
REM @echo off
REM 
REM echo. &
REM makeindex Thesis.glo -s Thesis.ist -t Thesis.glg -o Thesis.gls
REM echo. &
REM makeindex Thesis.nlo -s nomencl.ist -o Thesis.nls


REM -----------------------------------------------------------------
REM Full compilation
REM 
REM pdflatex Thesis
REM makeindex Thesis.glo -s Thesis.ist -t Thesis.glg -o Thesis.gls
REM makeindex Thesis.nlo -s nomencl.ist -o Thesis.nls
REM bibtex Thesis
REM pdflatex Thesis
REM pdflatex Thesis


REM -----------------------------------------------------------------
REM Full compilation (temporary files in subfolder \Build)
REM 
REM pdflatex -aux-directory=Build Thesis
REM makeindex Build\Thesis.glo -s Build\Thesis.ist -t Build\Thesis.glg -o Build\Thesis.gls
REM makeindex Build\Thesis.nlo -s Build\nomencl.ist -o Build\Thesis.nls
REM bibtex Build\Thesis
REM pdflatex -aux-directory=Build Thesis
REM pdflatex -aux-directory=Build Thesis


REM -----------------------------------------------------------------

