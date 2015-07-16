# Extract frequencies from given C2 word cloud data


# Read in all data --------------------------------------------------------
# Get character vector of all background terms
background <- read.csv('Data/Inflammation_C2_DistUp.csv')
background <- rownames(background)
background <- gsub('[A-Z]+_(.*)','\\1',background)
background <- unlist(strsplit(background,'_'))

# Words we do not want
C2_stopwords <- c("I", "HE", "THEIR", "ARE", "DOING", "YOU'VE", "WE'LL", "WOULDN'T", "WHEN'S", "UNTIL", "BEFORE", "UNDER", "BOTH", "SAME", "ME", "HIM", "THEIRS", "WAS", 
                  "WOULD", "WE'VE", "THEY'LL", "SHAN'T", "WHERE'S", "WHILE", "AFTER", "AGAIN", "EACH", "SO", "MY", "HIS", "THEMSELVES", "WERE", "SHOULD", "THEY'VE", 
                  "ISN'T", "SHOULDN'T", "WHY'S", "OF", "ABOVE", "FURTHER", "FEW", "THAN", "MYSELF", "HIMSELF", "WHAT", "BE", "COULD", "I'D", "AREN'T", "CAN'T", 
                  "HOW'S", "AT", "BELOW", "THEN", "MORE", "TOO", "WE", "SHE", "WHICH", "BEEN", "OUGHT", "YOU'D", "WASN'T", "CANNOT", "A", "BY", "TO", "ONCE", 
                  "MOST", "VERY", "OUR", "HER", "WHO", "BEING", "I'M", "HE'D", "WEREN'T", "COULDN'T", "AN", "FOR", "FROM", "HERE", "OTHER", "YOU", "IT", "THAT", 
                  "HAD", "SHE'S", "THEY'D", "HADN'T", "THAT'S", "BUT", "AGAINST", "IN", "WHERE", "NO", "YOUR", "ITS", "THESE", "HAVING", "IT'S", "I'LL", "DOESN'T", 
                  "WHO'S", "IF", "BETWEEN", "OUT", "WHY", "NOR", "OURSELVES", "HERSELF", "THIS", "HAS", "HE'S", "WE'D", "HAVEN'T", "LET'S", "AND", "ABOUT", "DOWN", 
                  "WHEN", "SUCH", "OURS", "HERS", "WHOM", "HAVE", "YOU'RE", "SHE'D", "HASN'T", "MUSTN'T", "THE", "WITH", "UP", "THERE", "SOME", "YOURS", "ITSELF", 
                  "THOSE", "DO", "WE'RE", "YOU'LL", "DON'T", "WHAT'S", "OR", "INTO", "ON", "HOW", "NOT", "YOURSELF", "THEY", "AM", "DOES", "THEY'RE", "HE'LL", 
                  "DIDN'T", "HERE'S", "BECAUSE", "THROUGH", "OFF", "ALL", "ONLY", "YOURSELVES", "THEM", "IS", "DID", "I'VE", "SHE'LL", "WON'T", "THERE'S", "AS", 
                  "DURING", "OVER", "ANY", "OWN", "LISTYOURUNIQUEWORDSHERE", "DN", "UP", "DN1", "UP1", "UP2", "DN2", "HR", "VS", "VIA", "EARLY", "LATE") 

background <- background[!background %in% C2_stopwords]

# Luminal A
luma <- read.csv('Data/C2_Inflammation_LumA_DistUp.txt',sep='\n')
luma <- as.character(luma[,1])
luma <- luma[!luma %in% C2_stopwords]

# Basal
basal <- read.csv('Data/C2_Inflammation_Basal_DistUp.txt',sep='\n')
basal <- as.character(basal[,1])
basal <- basal[!basal %in% C2_stopwords]

# Make a table with frequencies of all terms in each of the 3 data sets --------
terms <- background

bg.freq <- as.matrix(table(background),nrow=1)
bg.freq <- bg.freq/length(background)

luma.freq <- as.matrix(table(luma),nrow=1)
luma.freq <- luma.freq/length(luma)

basal.freq <- as.matrix(table(basal),nrow=1)
basal.freq <- basal.freq/length(basal)

common <- intersect(rownames(luma.freq),rownames(basal.freq))
bg.freq <- bg.freq[common,]
luma.freq <- luma.freq[common,]
basal.freq <- basal.freq[common,]

freq <- cbind(bg.freq,luma.freq,basal.freq)

# Direct Comparison between LumA and Basal --------------------------------
basal.tab <- table(basal)
basal.tab <- basal.tab[common]
luma.tab <- table(luma)
luma.tab <- luma.tab[common]
luma.vs.basal <- cbind(table(luma),table(basal))
colnames(luma.vs.basal) <- c('Luminal A','Basal')

pvals <- prop.test(luma.vs.basal)

# Adjustment for Background -----------------------------------------------

