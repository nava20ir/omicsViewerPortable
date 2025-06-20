SET ROPTS=--no-save --no-environ --no-init-file --no-restore --no-Rconsole
%cd%/App/R-Portable/bin/Rscript.exe %ROPTS% startApp.R 1> %cd%/run.log 2>&1
