# Average done with awk
BEGIN {
  num = 0;
  total = 0;
}
{
  if ($0 ~ /#.*/ || $0 ~ /@.*/) {
    print $0;
  } else {
    total = total + 1/$2;
    num = num + 1;
  }
}

END {
  print "Average is equal to " total/num;
}
