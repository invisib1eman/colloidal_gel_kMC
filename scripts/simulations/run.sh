salt_list=(0.9 1.7 3 9 80)
d_list=(0.919 0.638 0.459 0.275 0.0925)
step=0.1
N=8000
L=100
t=300
taucrawl=1
G=25
freeroll="True"
for i in "${!salt_list[@]}"; do
    echo "${salt_list[i]}"  # Correct reference to salt_list

    description="${salt_list[i]}mM_freeroll_step${step}_N_${N}_L_${L}_t_${t}_taucrawl_${taucrawl}"
    dump="../../trajectories/${description}.lammpsdump"
    log="../../logs/${description}.log"
    restart="../../restart/${description}.restart"
    data="../../data/${description}.data"

    colloidgel -N "${N}" -L "${L}" -G "${G}" -d "${d_list[i]}" --tcrawl "${taucrawl}" \
               -D "${description}" -m "${step}" -s "${t}" --dump "${dump}" \
               --log "${log}" --restart "${restart}" --data "${data}" &

done	
