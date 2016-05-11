args <- commandArgs(trailingOnly = TRUE)

require(seqinr,quietly=TRUE)  # install.packages("seqinr")
require(Matrix,quietly=TRUE)

# read first sequence in fasta file into character vector
readfasta = function(filename)
{
                                        # get string
    seq = read.fasta(filename,seqonly=TRUE, as.string=TRUE)[[1]]
                                        # get characters
    seq = strsplit(seq,split="")[[1]]
}

# convert character vector to vector of integers
tonumbers = function(letters)
{
    letters = as.numeric(factor(letters,levels=c('A','T','G','C')))
}

# show the first or second sequence, where state S corresonds to a '-' in this sequence
show.one = function(letters,path,S)
{
    row = c()
    x = 1
    for(i in 1:length(path))
    {
        if (path[i] == S)
        {
            row=c(row,'-')
        }
        else
        {
            row=c(row,letters[x])
            x = x + 1
        }
    }
    paste(row,collapse="")
}

# convert a sequence of states (path) into a readable alignment with letters
draw.a = function(letters1,letters2,path)
{
    paste(show.one(letters1,path,I),"\n",show.one(letters2,path,D),"\n",sep="",collapse="")
}


# wrap x on the range [0,b]
wrap0 = function(x, b)
{
    x = x %% (2*b)
    if (x > b) {
        x = 2*b - x
    }
    x
}

# wrap x on range [a,b]
wrap = function(x, a, b)
{
    a + wrap0(x-a,b-a)
}


# choose the 1st, 2nd, 3rd, etc. element of pr in proportion to their magnitude
choose = function(pr)
{
    pr = pr/sum(pr)
    u = runif(1)
    total = 0.0
    for(i in 1:length(pr))
    {
        total = total + pr[i]
        s = i
        if (u < total) {
            break
        }
    }
    s
}

# Choose integers for the states Match, Delete, Insert
M = 1
D = 2
I = 3
  

# Create the matrix of transition probabilities between states M,D,I
pairhmm = function(d,e)
{
    T = matrix(nrow=3,ncol=3)
    
    T[M,I] = d
    T[M,D] = d
    T[M,M] = 1-2*d
    
    T[I,M] = (1-e)*(1-2*d)
    T[I,I] = e + (1-e)*d
    T[I,D] = (1-e)*d
    
    T[D,M] = (1-e)*(1-2*d)
    T[D,I] = (1-e)*d
    T[D,D] = e + (1-e)*d

    T
}

# jukes-cantor rate matrix
jc = function()
{
}

# perform the forward algorithm to compute the dynamic programming matrix
forwardmatrix = function(seq1,seq2,T,P,pi)
{
}

# find the total probability of all alignments, from the forward matrix
total.prob = function(FM)
{
}

# perform backwards sampling to sample a random alignment in proportion to its probability
backsample = function(FM)
{
}

# read the input sequences
seq1letters = readfasta(args[1])
seq2letters = readfasta(args[2])
# how many iterations of MCMC
niter = args[3]

# convert the sequence strings to integers
seq1 = tonumbers(seq1letters)
seq2 = tonumbers(seq2letters)

Q = jc()
t = 0.25
P = expm(Q*t)
d=0.05
e=0.5

pi = rep(0.25, 4)
FM = forwardmatrix(seq1, seq2, pairhmm(d,e), P,pi)
path = backsample(FM)
write(draw.a(seq1letters,seq2letters,path), stderr())
