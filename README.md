Tree Swallow CMR model
================

# methods

## Study population

The study was conducted using a nest-box network established in 2004 in
southern Québec covering a 10 200 km2 area which is characterized by an
east-west gradient of agricultural intensification. We monitored 400
nest-boxes distributed equally among 40 farms spanning this gradient
every two days throughout the breeding season (April to August) from
2006 to 2019. Nests were monitored every two days during the breeding
season (from early May to mid-July). Females were captured in nest-boxes
during incubation and nestling rearing. At first capture, adult females
were marked individually with an official aluminium metal ring and
classified as second year (SY) or as after second year (ASY) according
to plumage coloration following Hussel (REF). Nestlings were marked at
12 days old and genetically sexed
\[@riouxpaquetteSevereRecentDecrease2014\]. \*\* include range of sample
size info \*\*

## Environment

We quantified local environment using several variables: mean
temperature (<sup>∘</sup>C), total precipitation, and days of cold snaps
during rearing period (defined as the period from 1st of june to 15th of
july), level of agricultural intensification and the level of
competition with it’s main competitor the house sparrow. Temperature was
measured using ibutton (REF) taking a recording at every hour with a
resolution of 0.5<sup>∘</sup>C placed on the underside of one nest boxes
on each farms. Precipitation was recorded every 2 days using
pluviometers located near a nest box on each farm and summed over the
rearing period. Cold snap was defined as the number of days during which
the maximum daily temperature measured by the ibutton did not exceed
18.8<sup>∘</sup>C, a measure which as previously been found to correlate
with fledging success in this population(Garret et al. 2021). To
quantify the level of agricultural intensification, we used a robust
principal component analysis based on the the compositional data.
Compositional data was assesed in situ each year within a radius of 500m
around each nest boxes and classified as forest, corn & soybean, forage
fields and cereals. We then calculated the mean percent cover of these
habitats across the 10 nest boxes on each farm and for each year
independently. The first component of the PCA explains 80.34% of the
variance in landcover, and has a strong positive correlation with corn
and soybean and negative with both forage fields and forest cover (see
Garret et al. 2021 for more details). House sparrow competition was
defined as the number of nestboxes in which at least 1 house sparrow
breeding attempt was made. All environmental variables were standardized
to a mean of 0 and sd of 1 prior to analysis.

## Vital rate estimation

Vital rates were estimated using a multi-state CMR model in order to
account for imperfect observation (Lebreton et al., 2009; Kéry & Schaub
2012). We considered 3 states: successfully reproductive,
non-successfully reproductive and dead ( reproductive success being
defined as having fledged at least 1 offspring). This allowed us to
control for the much higher capture probability of females with
successful reproduction than those without (Lagrange et al. XXXX). The
transitions between these states allowed us to estimate both survival
probability from one summer to the next and reproductive success. In
addition, the number of fledglings produced was jointly modeled using a
zero-inflated poisson model (REF). The model was implemented in nimble
(ref). MCMC chains were run for 60000 iterations with a burn-in of 50000
and a thinning of 5. Model convergence was assessed visually and using
Gelman and Rubin’s convergence diagnostic (Gelman & Rubin, 1992).

Table 1:

| parameter |                description                |
|:----------|:-----------------------------------------:|
| s         |           survival probability            |
| S         |                survival (0                |
| R         |             \> 1 fledgling (0             |
| r         | probability to have at least 1 fledgling  |
| F         |       Observed number of fledglings       |
| f         | conditional expected number of fledglings |
| i         |                individual                 |
| t         |                   year                    |
| l         |           location (farm, 1:40)           |
| A         |              effects of age               |
| a         |             age (N, SY, ASY)              |
| E         |      matrix of environment variables      |
| X         |    age-dependent effect of environment    |
| h         |   random effect of year (age\*response)   |
| j         |           random effect of farm           |
| k         |        random effect of individual        |

### survival and reproductive success

We estimated survival (*S* ∈ (0,1)) and reproductive success
(*R* ∈ (0,1)) using a hierarchical Bayesian multi-state model (Kéry and
Schaub 2012) with 3 states: 1) *R* = 0, 2) *R* = 1, 3) dead. Transitions
between state depended on both the probability to survival and the
probability to successfully produce fledgling (table 2). To simplify the
model, we assumed that there was no reproductive cost of reproduction
(*r* is the same when transiting from state 1 -> 2 and from 2 -> 2) and
no survival cost of reproduction (*s* is the same for state 1 and 2).
Both *s* and *r* were modeled as a function of age (3 age classes:
nestling, second year and after second year). Individual heterogeneity
in reproductive performance was accounted for using a random effects of
individual while un-explained annual variation in both reproductive
performance and survival was accounted for using year as a random
effect. Such that:

*l**o**g**i**t*(*s*<sub>*i*, *t*</sub>) = *A*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*s*</sup> + *E*<sub>*i*, *t*</sub> \* *X*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*s*</sup> + *h*<sub>*t*, *a*<sub>*i*, *t*</sub></sub> + *j*<sub>*l*<sub>*i*, *t*</sub>, 1</sub>

<!-- $$  -->
<!-- S_{i,t} \sim dbernoulli(s_{i,t}) -->
<!-- $$ -->

*l**o**g**i**t*(*r*<sub>*i*, *t*</sub>) = *A*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*r*</sup> + *E*<sub>*i*, *t* + 1</sub> \* *X*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*r*</sup> + *h*<sub>*t*, (*a*<sub>*i*, *t*</sub>+3)</sub> + *j*<sub>*l*<sub>*i*, *t* + 1</sub>, 2</sub> + *k*<sub>*i*, 1</sub>

<!-- $$ -->
<!-- R_{i,t} \sim dbernoulli(r_{i,t}) -->
<!-- $$ -->

Table 2: state-transition and observation matrix

|         |                        1: R=0                        |                      2: R=1                      |          3: dead           |
|:-------:|:----------------------------------------------------:|:------------------------------------------------:|:--------------------------:|
| 1: R=0  | *s*<sub>*i*, *t*</sub> \* (1−*r*<sub>*i*, *t*</sub>) | *s*<sub>*i*, *t*</sub> \* *r*<sub>*i*, *t*</sub> | 1 − *s*<sub>*i*, *t*</sub> |
| 2: R=1  | *s*<sub>*i*, *t*</sub> \* (1−*r*<sub>*i*, *t*</sub>) | *s*<sub>*i*, *t*</sub> \* *r*<sub>*i*, *t*</sub> | 1 − *s*<sub>*i*, *t*</sub> |
| 3: dead |                          0                           |                        0                         |             1              |
|  :—–:   |                         :—–:                         |                       :—–:                       |           :————:           |
|         |                      non-Breed                       |                      breed                       |          not seen          |
| 1: R=0  |                 *p*<sub>1, *t*</sub>                 |                        0                         |  1 − *p*<sub>1, *t*</sub>  |
| 2: R=1  |                          0                           |               *p*<sub>2, *t*</sub>               |  1 − *p*<sub>2, *t*</sub>  |
| 3: dead |                          0                           |                        0                         |             1              |

### Number of fledging

Number of fledgling follows a zero-inflated Poisson distribution.
Zero-inflation models can be considered as a mixture of 2 distribution
reflecting the probability to have more than 0 (state==2) and the number
of fledgling (REF). We can model this has:

*F*<sub>*i*, *t*</sub> ∼ *d**p**o**i**s**s**o**n*(*f*<sub>*i*, *t*</sub>\**R*<sub>*i*, *t*</sub>)

*R*<sub>*i*, *t*</sub> is already modeled by the multi-state CMR model,
so we only need to additionally model the *f*<sub>*i*, *t*</sub>. Just
like the reproductive success and survival, we model it as a function of
age, with id and (age\|yr) as random effects.

*e**x**p*(*f*<sub>*i*, *t* + 1</sub>) = *A*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*f*</sup> + *E*<sub>*i*, *t* + 1</sub> \* *X*<sub>*a*<sub>*i*, *t*</sub></sub><sup>*f*</sup> + *h*<sub>*t*, (*a*<sub>*i*, *t*</sub>+6)</sub> + *j*<sub>*l*<sub>*i*, *t*</sub>, 3</sub> + *k*<sub>*i*, 2</sub>

### random effects

### model selection

different candidate models were tested to reflect different effects of
environment.

| model |                  var                  |  k  | WAIC  | delta WAIC |
|:------|:-------------------------------------:|:---:|:-----:|:----------:|
| 1     |                  age                  |  9  | 10945 |    56.2    |
| 2     |             age\*coldSnap             | 18  | 10946 |    57.2    |
| 3     |           age\*(temp+prec)            | 27  | 10929 |    39.8    |
| 4     |               age\*hosp               | 18  | 10889 |     0      |
| 5     |                age\*AI                | 18  | 10914 |    25.4    |
| 6     |            age\*(AI+hosp)             | 27  | 10912 |    23.0    |
| 7     |         age\*(coldSnap+prec)          | 27  | 10943 |    54.4    |
| 8     |    Age \* (AI \* coldSnap \* prec)    | 63  | 10997 |   108.3    |
| 9     | Age \* (AI \* coldSnap \* prec +hosp) | 72  | 11027 |   137.8    |

### immigration rates

We estimate 3 immigration rates per year (see life cycle bellow): the
immigration of N to the SY class, of SY to ASY and of ASY to ASY. The
immigrating individual in the 2 later rates cannot be distinguished.
Based on some number from Esther about the proportion of returning
females which are SY, we assume that 48 % of arriving ASY where SY the
previous year (based on Carle-Pruneau et al. 2021). These numbers of
arriving individuals are then divided by the size of marked population
alive in the relevant stage the previous year (including those that are
alive, but not captured).

## demographic analysis

Leslie matrices were constructed using a postbreeding census approach
and assuming an offspring sex ratio of 0.5. ![](graphs/lifeCycle.png)

|     |                              N                               |                              SY                              |                             ASY                              |
|:----|:------------------------------------------------------------:|:------------------------------------------------------------:|:------------------------------------------------------------:|
| N   | *s*<sub>1</sub> \* *r*<sub>1</sub> \* *f*<sub>1</sub> \* 0.5 | *s*<sub>2</sub> \* *r*<sub>2</sub> \* *f*<sub>2</sub> \* 0.5 | *s*<sub>3</sub> \* *r*<sub>3</sub> \* *f*<sub>3</sub> \* 0.5 |
| SY  |              *s*<sub>1</sub> + *i*<sub>1</sub>               |                              0                               |                              0                               |
| ASY |                              0                               |              *s*<sub>2</sub> + *i*<sub>2</sub>               |              *s*<sub>3</sub> + *i*<sub>3</sub>               |

where *s*<sub>1</sub> is the survival from nestling to second year,
*r*<sub>1</sub> is the probability of a nestling to successfully fledge
a young next year (a second year) and *f*1 is the expected number of
fledgling for a successful second year; and so forth to the other age
classes. lambda sensitivities and elasticities were cacultated using the
popbio package in R. Sensitivities and elasticities of each fitness
component were decomposed from the demographic transition using the
approach by Van Tienderen (2000). These were calculated both for each
year of the study (by including the year random effect and mean annual
values for environmental variable) and at the population level (by
ignoring all random effects).

## variance partitioning of lambda

To quantify the relative importance of temporal fluctuation and spatial
variation in explaining the variation in fitness of the population, we
predicted fitness component for all age classes in a given environment,
build appropriate Leslie matrix and calculated lambda. Using a
bootstrapping approach,

# Results

based on WAIC selection, the best model only included the density of
house sparrow as explanatory variable.
