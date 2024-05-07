# Example cipher

- `cipher.sage` contains the implementation of the example cipher. It also exhaustively calculates the probability that one of the provided characteristics occur, given some predetermined round keys (can be used to check the output of `qdt.sage`).
- `qd-trails.sage` computes the quasidiferential transition matrix of the S-boxes, and uses path-finding to list all possible quasidifferential trails and their correlation for the given characteristics (also outputs the total).
