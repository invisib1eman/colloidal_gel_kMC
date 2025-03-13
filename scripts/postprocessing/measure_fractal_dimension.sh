salt_list=(0.9 1.7 3 9 80)
step=0.1
N=8000
L=100
t=300
# bin numbers
bins=30
#Unit length in A
R0=80
taucrawl_list=(1 100)
freeroll="True"
echo "salt step_size n_particles box_length total_time tau_crawl freeroll start_index end_index Df" >> ../../fractal_dimension/fractal_dimension.txt 
for salt in ${salt_list[@]}; do
	for taucrawl in ${taucrawl_list[@]}; do
	echo ${salt}
	python structural_fractal_dimension.py --salt ${salt} --step_size ${step} --n_particles ${N} --box_length ${L} --total_time ${t} --tau_crawl ${taucrawl} --freeroll ${freeroll} --bins ${bins} --UnitLength ${R0}
done
done
