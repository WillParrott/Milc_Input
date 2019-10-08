#!/bin/bash

# ran 840 on 16, 852 on 8 and 864 on 32 nodes


conf_start=864
conf_end=864



#conf_end=6294
dconf=12
#dconf=60

gen_input()
{
  local cfg=$1
  source /home/dc-parr2/venv3/bin/activate
  python3 ./in/WriteInput.py $cfg
}

for cfg in $(seq $conf_start $dconf $conf_end); do
  gen_input $cfg
  echo "configuration number: $cfg"
  sbatch ./submit/use_this_submit $cfg
done
