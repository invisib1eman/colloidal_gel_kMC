salt_list=(0.9 1.7 3 9 80)
step=0.1
N=8000
L=100
t=300
taucrawl_list=(1 100)
freeroll="True"
for salt in ${salt_list[@]}; do
	for taucrawl in ${taucrawl_list[@]}; do
	echo ${salt}
	python structural.py --salt ${salt} --step_size ${step} --n_particles ${N} --box_length ${L} --total_time ${t} --tau_crawl ${taucrawl} --freeroll ${freeroll}
done
done
