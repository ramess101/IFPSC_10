# Awk script to correct values found in viscosity output files
# argument 1: the correct box size; argument 2: the wrong box size

BEGIN { 
  factor=first / second;  # Correction factor
}
{
  if ($0 ~ /#.*/ || $0 ~ /@.*/) {
    print $0;
  } else {
    print $1 "   " $2*factor "  " $3*factor;
  }
}

