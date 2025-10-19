#!/bin/bash
# generate_combinations.sh

# Set values for each option
c_values=("1" "0.7")
d_values=("0" "1")
f_values=("default" "local")
l_values=("default" "50")
m_values=("high" "med" "low" "ultralow")
t_values=("4" "8")

# Basic command
echo "../FMAlign2x -i ../data/mt1x.fasta"

# -c
for c in "${c_values[@]}"; do
    echo "../FMAlign2x -i ../data/mt1x.fasta -c $c"
done

# -c -d
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d"
    done
done

# -c -d -f
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d -f $f"
        done
    done
done

# -c -d -f -l
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d -f $f -l $l"
            done
        done
    done
done

# -c -d -f -l -m
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for m in "${m_values[@]}"; do
                    echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d -f $f -l $l -m $m"
                done
            done
        done
    done
done

# -c -d -f -l -m -t
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for m in "${m_values[@]}"; do
                    for t in "${t_values[@]}"; do
                        echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d -f $f -l $l -m $m -t $t"
                    done
                done
            done
        done
    done
done

# all options
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for m in "${m_values[@]}"; do
                    for t in "${t_values[@]}"; do
                            echo "../FMAlign2x -i ../data/mt1x.fasta -c $c -d $d -f $f -l $l -m $m -t $t"
                    done
                done
            done
        done
    done
done
