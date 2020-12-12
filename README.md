# binpacking-hybrid-genetic-algorithm
A functional but incomplete hybrid genetic algorithm for the one-dimensional bin-packing problem.

There were several ambitions for this genetic algorithm that unfortunately were not sufficiently implemented within the timeframe for the project. Namely, there is a choice of "Lamarckian" and "Baldwinian" style of local search, but only Lamarckian is fully implemented. The backbone of the Baldwinian style is built, but the actual genetic algorithm was not coded to use it. The code is also not as clean as I'd like it to be.

The files that this works with are the bin-packing files from the JE Beasely's OR Library. The feature called "recycling" was used to work specifically with these files.

The code was built to support greedy or steep local search, but only greedy is used in the genetic algorithm code.

I used this project as an exploration of POSIX threads, and as such this will only compile on sufficiently POSIX-compatible systems. Importantly, there is a race condition when finding the best fitness and average fitness of the population. It does not impact the data or run of the algorithm and therefore I opted to not fix it, but it does result in inaccurate data recording.

**Features:**
- Permutation encoding of chromosomes
- First-fit heuristic decoder
- Order Crossover (OX)
- Tournament selection of size 2
- Full generational replacement except the best chromosome* from the previous generation which is kept
- Optional: "recycles" good permutations from previous problem instances
- Two types of local search: Random bin shuffling and random element swaps
  - Random bin shuffling swaps two bins in a decoded solution, then re-encodes it back into a first-fit-compatible permutation, then decodes it once again
  - Random element swaps simply swaps two elements in the permutation and decodes it
- Optional: Use a hill-climbing procedure to generate initial population
- Optional: Use local search in place of mutation
- Prints out average and best fitnesses** alongside the amount of time spent in each generation

\*- Race condition causes a really good, but not the best, chromosome to be kept sometimes

\*\*- Race condition makes the average and best data to be somewhat inaccurate
