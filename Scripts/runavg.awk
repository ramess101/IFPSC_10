# Running Average done with awk
BEGIN {
  num = 0;
  total = 0;
}

{
  if (NR > 100000) {  # Cut the junk off the beginning
    if ($0 ~ /#.*/ || $0 ~ /@.*/) {
      print $0;
    } else {
      total = total + 1/$2;
      num = num + 1;
      print $1 "	" total/num
    }
  } else {
    print $1 "	" 0
  }
}

END {
  print total/num > "visco_avg.txt"
}
