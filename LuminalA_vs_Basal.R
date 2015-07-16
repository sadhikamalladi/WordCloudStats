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

# Frequency table
luma.freq <- unlist(table(luma))/length(luma)
basal.freq <- unlist(table(basal))/length(basal)
bg.freq <- unlist(table(background))/length(background)

freq <- matrix(NA,nrow=length(bg.freq),ncol=3)
freq[,1] <- bg.freq
rownames(freq) <- names(bg.freq)
freq[names(luma.freq),2] <- luma.freq
freq[names(basal.freq),3] <- basal.freq
colnames(freq) <- c('Background','Luminal A','Basal')

# Counts table
luma.ct <- unlist(table(luma))
basal.ct <- unlist(table(basal))
bg.ct <- unlist(table(background))

count <- matrix(NA,nrow=length(bg.ct),ncol=3)
count[,1] <- bg.ct
rownames(count) <- names(bg.ct)
count[names(luma.ct),2] <- luma.ct
count[names(basal.ct),3] <- basal.ct
colnames(count) <- c('Background','Luminal A','Basal')

# prop test
pvals <- rep(1,nrow(count))
names(pvals) <- rownames(count)
for (term in rownames(count)) {
  luma.yes <- count[term,2]
  basal.yes <- count[term,3]
  luma.no <- length(luma) - luma.yes
  basal.no <- length(basal) - basal.yes
  
  mat <- matrix(NA, nrow=2, ncol=2)
  mat[1,1] <- luma.yes
  mat[1,2] <- luma.no
  mat[2,1] <- basal.yes
  mat[2,2] <- basal.no
  
  p <- prop.test(mat)$p.value
  pvals[term] <- p
}


