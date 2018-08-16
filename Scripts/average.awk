# Awk average
BEGIN {
  num = 0;
  total = 0;
}
{
  total = total + $1
  num = num + 1
}
END {
  print "Average is equal to " total /num;
}
