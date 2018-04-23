#!/bin/csh

setenv HOME /phenix/u/$LOGNAME
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.cshrc
source $HOME/.login

./test_prod.exe $1

