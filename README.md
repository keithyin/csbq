* default model: model/default.txt

```
# for hsc 
csbq --mode hsc -q reads.bam/fa/fq -t contig.fa/fq -o cali.fq

# for smc
csbq --mode smc -q sbr.bam -t smc.bam -o cali.bam

```


## Error assignment

### Mismatch

```
           *
    smc: ACGTATA
    sbr: ACATATA
```


### Deletion

```
           *
    smc: ACATATA
    sbr: AC-TATA
  
```


### Insertion

```
            *
    smc: AC-TATA
    sbr: ACGTATA
  
```


### Homopolymer Insertion or Deletion

Indels should always be left-aligned and the error probability is only given for the first base in a homopolymer.

```
           *                  *
   smc: AC-GGGTATA     smc: ACGGGGTATA 
   sbr: ACGGGGTATA     sbr: AC-GGGTATA
   
```


# UnigramModel

If the model path is not specified, the default model parameters are as follows

* AC 0.9 means that, given the ref base is A, the probability of calling C is 0.9
* A 0.25 means that, in the current reference genome, the probability of A appearing is 0.25.

```
AA 0.9
AC 0.02
AG 0.02
AT 0.02
A- 0.02
A+ 0.02
CA 0.02
CC 0.9
CG 0.02
CT 0.02
C- 0.02
C+ 0.02
GA 0.02
GC 0.02
GG 0.9
GT 0.02
G- 0.02
G+ 0.02
TA 0.02
TC 0.02
TG 0.02
TT 0.9
T- 0.02
T+ 0.02
A 0.25
C 0.25
G 0.25
T 0.25
```
