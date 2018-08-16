# Awk script
BEGIN { print "# Regardless of what appears below, this file contains time and NEMD viscosity" }
{
  if ($0 ~ /#.*/ || $0 ~ /@.*/) {
    print $0;
  } else {
    if ($2 == 0) {
      print 0
    } else {
      print $1 "   " 1.0/$2;
    }
  }
}
