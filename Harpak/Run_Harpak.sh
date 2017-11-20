#!/bin/bash
sbatch -p long --exclusive -o Harpak-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Harpak_Guess_1.sh  
sbatch -p long --exclusive -o Harpak-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Harpak_Guess_2.sh  
sbatch -p long --exclusive -o Harpak-%j.out --mail-type=FAIL --mail-user=xji3@ncsu.edu ./ShFiles/Harpak_Guess_3.sh  
