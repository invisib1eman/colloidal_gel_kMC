mamba activate porous
salt_list=(1.7 3 9 80)
step=0.1
N=8000
L=100
t=300
taucrawl=100
freeroll="True"
for salt in ${salt_list[@]}; do
	echo ${salt}
	python structural.py --salt ${salt} --step_size ${step} --n_particles ${N} --box_length ${L} --total_time ${t} --tau_crawl ${taucrawl} --freeroll ${freeroll}
done
