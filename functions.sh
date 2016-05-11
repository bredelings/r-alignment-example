function Rexec {
  local script=$1;
  shift;
  R --slave --vanilla --args "$@" < $script
}
