#!/bin/bash
# generate_combinations.sh

# Definir valores para cada opção
c_values=("1" "0.7")
d_values=("0" "1")
f_values=("default" "local")
l_values=("default" "50")
p_values=("mafft" "halign3")
t_values=("1" "2")
x_values=("0" "1")

# Comando básico
echo "./FMAlign2x -i data/mt1x.fasta"

# Combinações com -c
for c in "${c_values[@]}"; do
    echo "./FMAlign2x -i data/mt1x.fasta -c $c"
done

# Combinações com -c -d
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d"
    done
done

# Combinações com -c -d -f
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d -f $f"
        done
    done
done

# Combinações com -c -d -f -l
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d -f $f -l $l"
            done
        done
    done
done

# Combinações com -c -d -f -l -p
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for p in "${p_values[@]}"; do
                    echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d -f $f -l $l -p $p"
                done
            done
        done
    done
done

# Combinações com -c -d -f -l -p -t
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for p in "${p_values[@]}"; do
                    for t in "${t_values[@]}"; do
                        echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d -f $f -l $l -p $p -t $t"
                    done
                done
            done
        done
    done
done

# Combinações completas com todas as opções
for c in "${c_values[@]}"; do
    for d in "${d_values[@]}"; do
        for f in "${f_values[@]}"; do
            for l in "${l_values[@]}"; do
                for p in "${p_values[@]}"; do
                    for t in "${t_values[@]}"; do
                        for x in "${x_values[@]}"; do
                            echo "./FMAlign2x -i data/mt1x.fasta -c $c -d $d -f $f -l $l -p $p -t $t -x $x"
                        done
                    done
                done
            done
        done
    done
done
