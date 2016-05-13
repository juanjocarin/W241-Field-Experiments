## @knitr Load_Libraries
library(boot)
library(car)
library(complmrob)
library(dgof)
library(effsize)
library(ggplot2)
library(gmodels)
library(grid)
library(gridExtra)
library(kfigr)
library(knitr)
library(MASS)
library(memisc)
library(nlme)
library(pander)
library(plyr)
library(png)
library(pwr)
library(reshape)
library(robust)
library(robustbase)
library(simpleboot)
library(xtable)




## @knitr Definition_Functions
sig_stars = function(p) {
	stars = symnum(p, na = F, cutpoints = c(0, .001, .01, .05, .1, 1), 
						symbols=c("**`***`**","**`** `**", "**`*  `**", "**.  **", "   "))
	return(stars)
}


# A function to calculate RSEs
RSEs <- function(model){
	require(sandwich, quietly = TRUE)
	require(lmtest, quietly = TRUE)
	newSE <- vcovHC(model)
	coeftest(model, newSE)
}


# Function to calculate CSEs
cl <- function(fm, cluster){
	require(sandwich, quietly = FALSE)
	require(lmtest, quietly = FALSE)
	M <- length(unique(cluster))
	N <- length(cluster)
	K <- fm$rank
	dfc <- (M/(M-1))*((N-1)/(N-K))
	uj <- apply(estfun(fm),2, function(x) tapply(x, cluster, sum));
	vcovCL <- dfc*sandwich(fm, meat=crossprod(uj)/N)
	coeftest(fm, vcovCL)
}


# A function to detect outliers based on the median
find_outliers_mdn <- function(vector, lower, upper) {
	qtile <- as.numeric(quantile(vector, probs = c(lower/100, 1 - upper/100)))
	outliers <- vector < qtile[1] | vector > qtile[2]	
	return(vector[outliers])
}


# A function to detect outliers based on the mean
find_outliers_mn <- function(vector, lower, upper) {
	outliers <- scale(vector) < -qnorm(1 - lower/100) | 
		scale(vector) > qnorm(1 - upper/100)
	return(vector[outliers])
}


# A function to plot mosaics
ggMMplot <- function(var1, var2){
	levVar1 <- length(levels(var1))
	levVar2 <- length(levels(var2))
	jointTable <- prop.table(table(var1, var2))
	plotData <- as.data.frame(jointTable)
	plotData$marginVar1 <- prop.table(table(var1))
	plotData$var2Height <- plotData$Freq / plotData$marginVar1
	plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) + 
		plotData$marginVar1 / 2
	ggplot(plotData, aes(var1Center, var2Height)) + 
		geom_bar(stat = "identity", aes(width = marginVar1, fill = var2), 
					colour = "Black") + 
		geom_text(aes(label = as.character(var1), x = var1Center, y = 1.05)) + 
		xlab("Group")
}


# A function to perform regression based on bootstrapped data
bootReg <- function(formula, data, i) {
	d <- data[i, ]
	fit <- lm(formula, data = d)
	return(coef(fit))
}


# A function to estimate p-values (based on robust SEs)
estimate_pvalues <- function(df) {
	pvalues <- RSEs(lm(Score ~ Group, df))[c(2:3), 4]
	names(pvalues) <- substring(names(pvalues), 6)
	return(pvalues)
}


# A function to simulate a dataset with same means and SDs
simulate_study <- function(df, N) {
	mean_Score <- as.numeric((with(df, by(df, Group, function(x) 
		mean(x$Score)))))
	sd_Score <- as.numeric((with(df, by(df, Group, function(x) 
		sd(x$Score)))))
	sizes <- round(as.numeric((with(df, by(df, Group, nrow)))) / nrow(df) * N)
	df_sim <- data.frame(Group = c(rep("Control", sizes[1]), 
											 rep("Loser", sizes[2]), 
											 rep("Winner", sizes[3])), 
								Score = c(rnorm(sizes[1], mean_Score[1], sd_Score[1]), 
											 rnorm(sizes[2], mean_Score[2], sd_Score[2]), 
											 rnorm(sizes[3], mean_Score[3], sd_Score[3])))
	return(df_sim)
}


# A function to estimate statistical power
estimate_power <- function(df, N) {
	results <- replicate(1e3, estimate_pvalues(simulate_study(df, N)))
	results <- as.data.frame(t(results))
	effect.detected <- results[, 1] < 0.05/2 | results[, 2] < 0.05/2
	return(mean(effect.detected))
}


calculate_mean <- function(df, variable1, variable2) {
	formula <- ifelse(!missing(variable2), 
							paste0(".d$", variable1, " - ", ".d$", variable2, " ~ 1"), 
							paste0(".d$", variable1, " ~ 1"))
	mean_score <- t(ddply(df, .(Group), function(.d) 
		data.frame(x = paste0(formatC(cl(lm(as.formula(formula)), 
													.d$Cluster_ID)[, 1], digits = 2, 
												format = "f", drop0trailing = FALSE), 
									 " (", formatC(cl(lm(as.formula(formula)), 
									 					  .d$Cluster_ID)[, 2], digits = 2, 
									 				  format = "f", 
									 				  drop0trailing = FALSE), ")"))))
	total <- ddply(df, .(), function(.d) 
		data.frame(x = paste0(formatC(cl(lm(as.formula(formula)), 
													.d$Cluster_ID)[, 1], digits = 2, 
												format = "f", drop0trailing = FALSE), 
									 " (", formatC(cl(lm(as.formula(formula)), 
									 					  .d$Cluster_ID)[, 2], digits = 2, 
									 				  format = "f", drop0trailing = FALSE), 
									 ")")))
	means <- c(mean_score[2, ], as.character(total[1, 2]))
	names(means) <- c(mean_score[1, ], "TOTAL")
	return(means)
}


plot_density <- function(df) {
	Mean_Score_Group <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score)))
	colnames(Mean_Score_Group)[2] <- "Score"
	CIlower <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score)-1.96*sd(.d$Score)/sqrt(length(.d$Score))))
	colnames(CIlower)[2] <- "CIlower"
	CIupper <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score)+1.96*sd(.d$Score)/sqrt(length(.d$Score))))
	colnames(CIupper)[2] <- "CIupper"
	density_Score <- ggplot(df, aes(Score))
	density_Score + geom_density(fill = "White") + 
		ylab("Frequency") + xlab("Score in the 2nd test") + 
		geom_vline(data = Mean_Score_Group, 
					  aes(xintercept = Score, colour = Group), size = 1) + 
		geom_vline(data = CIlower, aes(xintercept = CIlower, colour = Group), 
					  size = 0.5, linetype = "dashed") + 
		geom_vline(data = CIupper, aes(xintercept = CIupper, colour = Group), 
					  size = 0.5, linetype = "dashed") + 
		facet_grid(Group ~ .)
}


plot_density_DIFF <- function(df) {
	Mean_Score_DIFF_Group <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score - .d$Score_PRE)))
	colnames(Mean_Score_DIFF_Group)[2] <- "Score_DIFF"
	CIlower <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score - .d$Score_PRE) - 
					  	1.96*sd(.d$Score - .d$Score_PRE)/sqrt(length(.d$Score))))
	colnames(CIlower)[2] <- "CIlower"
	CIupper <- ddply(df, .(Group), function(.d) 
		data.frame(x = mean(.d$Score - .d$Score_PRE) + 
					  	1.96*sd(.d$Score - .d$Score_PRE)/sqrt(length(.d$Score))))
	colnames(CIupper)[2] <- "CIupper"
	density_Score <- ggplot(df, aes(Score - Score_PRE))
	density_Score + geom_density(fill = "White") + 
		ylab("Frequency") + xlab("Score difference") + 
		geom_vline(data = Mean_Score_DIFF_Group, 
					  aes(xintercept = Score_DIFF, colour = Group), size = 1) + 
		geom_vline(data = CIlower, aes(xintercept = CIlower, colour = Group), 
					  size = 0.5, linetype = "dashed") + 
		geom_vline(data = CIupper, aes(xintercept = CIupper, colour = Group), 
					  size = 0.5, linetype = "dashed") + 
		facet_grid(Group ~ .)
}


plot_line <- function(df) {
	df_2 <- df[, c("User_ID", "Gender", "Score_PRE", "Score", "Group")]
	df_2 <- melt(df_2, id = c("User_ID", "Gender", "Group"), 
					 measured = c("Score_PRE", "Score"))
	line_Score <- ggplot(df_2, aes(variable, value, colour = Group), 
								environment = environment())
	line_Score + stat_summary(fun.y = mean, geom = "point") + 
		stat_summary(fun.y = mean, geom = "line", aes(group = Group)) + 
		stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + 
		ylab("Score") + scale_x_discrete("", labels = c("1st test", "2nd test"))
}




## @knitr Open_Dataset
Results <- read.csv("final_results.csv")
nrow(Results)
head(Results)
row.names(Results) <- NULL




## @knitr Clean_Dataset_1
# First we check that there are no missing values regarding to ID, gender, age, 
	#education level, math confidence and IP address.
all(!is.na(Results[, 1:6]))
# TRUE, so the 363 records contain data about those 6 variables.

# Then we look for different user IDs that may correspond to the same 
	# participant.
unique_fields <- c("IP_Address", "Gender", "Age")

# Results_vector <- apply(Results[, unique_fields], 1 , function(x) 
# 	paste(x, collapse="_"))
# Results_table <- table(Results_vector)
# same_participant <- Results[which(Results_vector %in% 
# 											 	names(which(Results_table > 1))), c(1:6)]
# same_participant <- same_participant[order(same_participant$IP_Address, 
# 														 same_participant$Age, 
# 														 same_participant$User_ID), ]

same_participant <- Results[duplicated(Results[, unique_fields]) | 
									 	duplicated(Results[, unique_fields], 
									 				  fromLast = TRUE), c(1:6)]
same_participant <- same_participant[order(same_participant$IP_Address, 
														 same_participant$Age, 
														 same_participant$User_ID), ]
same_participant_IDs <- unique(lapply(same_participant$IP_Address, function(x) 
	same_participant$User_ID[same_participant$IP_Address == x]))
nrow(same_participant)
length(same_participant_IDs)
# There are 39 IDs that may correspond to only 17 participants.

# But...
# ID 120 has different math confidence (2) than IDs 127 and 131 (1)
# ID 117 has different education level (7) than ID 178 (6)
# Maybe they are really different individuals.
# Let's check who have same values of gender, age, and IP address, but different 
	# education level or math confidence.
# We could have searched for duplicates of the 5 variables by including the last 
	# 2 variables in "unique_fields."
diff_Edu_Math <- logical(0)
for (i in (1:length(same_participant_IDs))) {
	diff_Edu_Math <- c(diff_Edu_Math, 
							 !all(duplicated(same_participant$Math_Conf[
							 	same_participant$User_ID %in% 
							 		same_participant_IDs[[i]]])[-1]) | 
							 	!all(duplicated(same_participant$Edu_Level[
							 		same_participant$User_ID %in% 
							 			same_participant_IDs[[i]]])[-1]))
}
length(which(diff_Edu_Math == TRUE))
same_participant[same_participant$User_ID %in% 
					  	unlist(same_participant_IDs[which(diff_Edu_Math == TRUE)]), ]
# Education level or math confidence differ in 5 of the 17 subsets of 
	# duplicates... but only by 1 unit, so we will assume the participant entered 
	# a wrong value the first time he or she did it.

# We consider we have enough data from a participant (and hence that he or she 
	# has finished the experiment) if we have collected at least the main three 
	# variables: pre-test score, group, and score. If any of them is missing. We 
	# consider the participant attrited.

# For those subjects who participated more than once, we will only admit the 
	# first time they finished the experiment (any previous---unfinished---
	# attempt, as well as any subsequent attempt---whether they finished or 
	# not---is discarded).
# If one of those participants never finished the experiment, we consider his or 
	# her last attempt as attrition.
needed_fields <- c("Score_PRE", "Group", "Score")
Valid_IDs <- Non_Valid_IDs <- rep(NA, length(same_participant_IDs))
for (i in (1:length(same_participant_IDs))) {
	for (j in same_participant_IDs[[i]]) {
		if (all(!is.na(Results[Results$User_ID == j, needed_fields])) & 
			 	is.na(Valid_IDs[i])) {Valid_IDs[i] <- j}
	}
	if (is.na(Valid_IDs[i])) {Non_Valid_IDs[i] <- j}
}
Valid_IDs <- Valid_IDs[!is.na(Valid_IDs)]
Non_Valid_IDs <- Non_Valid_IDs[!is.na(Non_Valid_IDs)]
same_participant_IDs <- unlist(same_participant_IDs)
Discarded_IDs <- same_participant_IDs[-which(same_participant_IDs %in% 
																c(Valid_IDs, Non_Valid_IDs))]
length(Discarded_IDs)

Results <- Results[-which(Results$User_ID %in% Discarded_IDs), ]
Results0 <- Results
nrow(Results)




## @knitr Check_Randomization_1
# Before continuing with the analysis, we will check whether the random 
	# allocation of experimental conditions was done properly.
# We don't have to worry about the duplicates we have omitted, since the website 
	# already assigned the same group to all IDs sharing IP address, age and 
	# gender
Group <- factor(Results0$Group, levels = c(0:2), 
					 labels = c("Control", "Loser", "Winner"))
CrossTable(Group)

prop <- CrossTable(Group)
prop <- matrix(c(prop$t, prop$prop.row), ncol = 2)
prop <- rbind(prop, apply(prop, 2, sum))
colnames(prop) <- c("No. Participants", "Proportion")
rownames(prop) <- c(levels(Group), "TOTAL")
prop[, 2] <- paste0(formatC(100*prop[, 2], digits = 2, format = "f", 
									 drop0trailing = FALSE), "%")
# 37% (rather than 33.3%) of the subjects were assigned to the Losers' group.
# So let's check if the distribution is uniform (p = 1/3) or not:
#Kolmogorov-Smirnov Test (for discrete distributions!!)
(KS <- dgof::ks.test(Results0$Group, ecdf(0:2)))
dgof::ks.test(Results0$Group, ecdf(0:2), simulate.p.value = TRUE, B=1e3)
# Chi-Square Test
chi2 <- chisq.test(x = table(Results$Group), p = rep(1/3, 3))
chisq.test(x = table(Results$Group), p = rep(1/3, 3), simulate.p.value = TRUE, 
			  B = 1e3)




## @knitr Check_Randomization_2
kable(prop, align = "r", caption = "Subject Distribution by Group")

# qplot(Results0$Group, stat = "ecdf", geom = "step", size = I(1)) + 
# 	scale_y_continuous("F(x)", breaks = c(0,1/3,2/3,1), 
# 							 labels = c("0", "1/3", "2/3", "1")) + 
# 	scale_x_continuous( "Group", breaks = c(0:2), 
# 							  labels = levels(Group)) +
# 	theme(panel.grid.major.y = element_line(colour = "orange", linetype = 2))




## @knitr Clean_Dataset_2
# Now that we have dealt with the duplicates, we focus on invalid (very young) 
	# participants
Results <- Results[Results$Age >= 18, ]
nrow(Results)
# This is our cleaned dataset, that has 339 records
row.names(Results) <- NULL
write.csv(Results, "Results_cleaned.csv", na = "", row.names = FALSE)

# Next we convert categorical variables to factors
Results$Group <- factor(Results$Group, levels = c(0:2), 
								labels = c("Control", "Loser", "Winner"))
# For Education and Math Confidence levels, the baseline will be the most 
	# common one (we do this because the lowest level may be so rare that is 
	# not present in some Group)
Results$Edu_Level <- as.factor(Results$Edu_Level)
most_common_EduLevel <- which(table(Results$Edu_Level) == 
											max(table(Results$Edu_Level)))
if (most_common_EduLevel != 1) {
	Results$Edu_Level <- factor(Results$Edu_Level, 
										 levels = levels(Results$Edu_Level)[
										 	c(most_common_EduLevel:9, 
										 	  (most_common_EduLevel-1):1)])
}
Results$Math_Conf <- as.factor(Results$Math_Conf)
most_common_MathConf <- which(table(Results$Math_Conf) == 
											max(table(Results$Math_Conf)))
if (most_common_MathConf != 1) {
	Results$Math_Conf <- factor(Results$Math_Conf, 
										 levels = levels(Results$Math_Conf)[
										 	c(most_common_MathConf:9, 
										 	  (most_common_MathConf-1):1)])
}
# We will also use the most common gender as baseline:
Results$Gender <- factor(Results$Gender, 
								 levels = levels(Results$Gender)[c(2:1)])




## @knitr Check_Balanced_Groups_1
# Given random assignment of experimental condition (Group), the covariates 
	# should be well balanced across Groups.
df <- Results

cl(lm(as.numeric(Gender) - 1  ~ Group, df), df$Cluster_ID)
CrossTable_Gender <- CrossTable(df$Gender, df$Group, fisher = TRUE, 
										  expected = TRUE, prop.chisq = FALSE, 
										  prop.r = FALSE, prop.t = FALSE)
(balance_Gender <- CrossTable(df$Gender, df$Group, fisher = TRUE, 
										expected = TRUE, prop.chisq = FALSE, 
										prop.r = FALSE, prop.t = FALSE))
kable(round(100*as.matrix(ftable(prop.table(table(df$Gender, df$Group), 2))), 
				1), align = "r")
# mosaicplot(xtabs(~ Group + Gender, df))
# ggMMplot(df$Group, df$Gender) + ylab("Gender") + 
# 	scale_fill_discrete(name="Gender")
hist_Gender <- ggplot(df, aes(Group, fill = Gender))
y0 <- cumsum(prop.table(table(df$Gender))); names(y0) <- NULL
hist_Gender + geom_bar(position="fill") + ylab("Proportion") + 
	ggtitle("Proportion of genders in each Group") + 
	scale_fill_discrete(name = "Gender", guide = guide_legend(reverse=TRUE)) + 
	geom_hline(aes(yintercept = y0[1]), colour = "Black", linetype = "dashed", 
				  size = 0.5)
# Gender is equally distributed among the 3 experimental conditions (there are 
	# no significant differences between the percentage of men and women in each
	# Group)

cl(lm(Age ~ Group, df), df$Cluster_ID)
boxplot_Age <- ggplot(df, aes(factor(0), Age))
boxplot_Age + geom_boxplot() +  ylab("Age") + 
	scale_x_discrete("") +   
	ggtitle("Boxplot of Ages in each Group") + facet_grid(. ~ Group) + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
hist_Age <- ggplot(df, aes(Age))
hist_Age + geom_histogram(binwidth = 5, fill = "White", colour = "Black", 
								  aes(y = ..density..)) + ylab("Frequency") + 
	xlab("Age") + ggtitle("Histogram of participants' age in each group") + 
	facet_grid(Group ~ .)
# df$Age_Group <- as.factor(recode(df$Age, "18:29 = '18..29'; 
# 												  30:45 = '30..45'; 46:65 = '45..65'; 
# 												  66:99 = '>65'"))
# df$Age_Group <- factor(df$Age_Group, 
# 									 levels = levels(df$Age_Group)[c(2:4, 1)])
Q <- quantile(df$Age)
df$Age_Group <- as.factor(car::recode(df$Age, 
											paste0(Q[1], ":", Q[2], " = '", 
													 Q[1], "..", Q[2], "'; ", Q[2]+1, ":", 
													 Q[3], " = '", Q[2]+1, "..", Q[3], 
													 "'; ", Q[3]+1, ":", Q[4], " = '", 
													 Q[3]+1, "..", Q[4], "'; ", Q[4]+1, 
													 ":", Q[5], " = '", Q[4]+1, "..", Q[5], 
													 "'")))
cl(lm(as.numeric(Age_Group) ~ Group, df), df$Cluster_ID)
Age_Proportion <- lapply(levels(df$Age_Group), function(x) 
	cl(lm(Age_Group == x ~ Group, df), df$Cluster_ID))
names(Age_Proportion) <- levels(df$Age_Group)
Age_Proportion
CrossTable(df$Age_Group, df$Group, expected = TRUE, prop.chisq = FALSE, 
			  prop.r = FALSE, prop.t = FALSE)
kable(round(100*as.matrix(prop.table(table(df$Age_Group, df$Group), 2)), 1), 
		align = "r")
hist_AgeGroup <- ggplot(df, aes(Group, fill = Age_Group))
y0 <- cumsum(prop.table(table(df$Age_Group))); names(y0) <- NULL
hist_AgeGroup + geom_bar(position="fill") + ylab("Proportion") + 
	ggtitle("Proportion of age groups in each Group") + 
	scale_fill_discrete(name = "Age group", guide = guide_legend(reverse=TRUE)) + 
	geom_hline(aes(yintercept = y0[1]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[2]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[3]), colour = "Black", linetype = "dashed", 
				  size = 0.5)
# Age does not predict which Group a subject was assigned, either.

df$Edu_Level2 <- factor(df$Edu_Level, levels = c(2:9))
cl(lm(as.numeric(Edu_Level2)+1 ~ Group, df), df$Cluster_ID)
EduLevel_Proportion <- lapply(levels(df$Edu_Level), function(x) 
	cl(lm(Edu_Level == x ~ Group, df), df$Cluster_ID))
names(EduLevel_Proportion) <- levels(df$Edu_Level)
EduLevel_Proportion
CrossTable(df$Edu_Level, df$Group, expected = TRUE, prop.chisq = FALSE, 
			  prop.r = FALSE, prop.t = FALSE)
kable(round(100*as.matrix(prop.table(table(df$Edu_Level, df$Group), 2)), 1), 
		align = "r")
hist_EduLevel <- ggplot(df, aes(Group, fill = Edu_Level))
y0 <- cumsum(prop.table(table(df$Edu_Level))); names(y0) <- NULL
hist_EduLevel + geom_bar(position="fill") + ylab("Proportion") + 
	ggtitle("Proportion of education levels in each Group") + 
	scale_fill_discrete(name="Education\nlevel", 
							  guide = guide_legend(reverse=TRUE)) + 
	geom_hline(aes(yintercept = y0[1]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[2]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[3]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[4]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[5]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[6]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[7]), colour = "Black", linetype = "dashed", 
				  size = 0.5)
# There are some (significant) differences abmong Groups in the education level
	# (e.g., there are few subjects assigned to the Losers' Group with an
	# education level of 5---attended college but not finished studies---; but 
	# many education levels were so infrequent in our sample that differences
	# across Groups can be considered as normal).
# But overall, group does not seem to predict the education level.

cl(lm(as.numeric(Math_Conf) ~ Group, df), df$Cluster_ID)
MathConf_Proportion <- lapply(levels(df$Math_Conf), function(x) 
	cl(aov(Math_Conf == x ~ Group, df), df$Cluster_ID))
names(MathConf_Proportion) <- levels(df$Math_Conf)
MathConf_Proportion
CrossTable(df$Math_Conf, df$Group, fisher = TRUE, expected = TRUE, 
			  prop.chisq = FALSE, prop.r = FALSE, prop.t = FALSE)
kable(round(100*as.matrix(prop.table(table(df$Math_Conf, df$Group), 2)), 2), 
		align = "r", row.names = TRUE)
hist_MathConf <- ggplot(df, aes(Group, fill = Math_Conf))
y0 <- cumsum(prop.table(table(df$Math_Conf))); names(y0) <- NULL
hist_MathConf + geom_bar(position="fill") + ylab("Proportion") + 
	ggtitle("Proportion of math confidence levels in each Group") + 
	scale_fill_discrete(name="Math\nconfidence", 
							  guide = guide_legend(reverse=TRUE)) + 
	geom_hline(aes(yintercept = y0[1]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[2]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[3]), colour = "Black", linetype = "dashed", 
				  size = 0.5) + 
	geom_hline(aes(yintercept = y0[4]), colour = "Black", linetype = "dashed", 
				  size = 0.5)
# There is not any significant difference among Groups in the math confidence 
	# level (e.g., there are many subjects assigned to the Losers' Group with a
	# math confidence level of 1---they find arithmetic operations easy---; 
	# subjects with math confidence levels from 3 to 5 are so rare that 
	# differences across Groups can be considered as normal).


## @knitr Check_Balanced_Groups_2
create_regtable <- function(model, df) {
	model2 <- summary(model)
	model3 <- cl(model, df$Cluster_ID)
	estimate <- unlist(lapply(c(2:3, 1), function(x) 
		paste0(formatC(model3[x, 1], digits = 3, format = "f", 
							drop0trailing = FALSE), sig_stars(model3[x, 4]))))
	SE <- unlist(lapply(c(2:3, 1), function(x) 
		paste0("(", formatC(model3[x, 2], digits = 3, format = "f", 
								  drop0trailing = FALSE), ")  ")))
	N <- paste0(nrow(df), "   ")
	R2 <- paste0(formatC(model2$r.squared, digits = 3, format = "f", 
								drop0trailing = FALSE), "   ")
	Fstatistic <- paste0(formatC(model2$fstatistic[1], digits = 3, format = "f", 
							  drop0trailing = FALSE), "   ")
	pvalue <- paste0(formatC(1 - pf(model2$fstatistic[1], 2, 300), digits = 3, 
									 format = "f", drop0trailing = FALSE), "   ")
	table <- matrix(c(t(matrix(c(estimate, SE), ncol = 2)), R2, Fstatistic, 
							pvalue, N), ncol = 1)
	rownames(table) <- c("**Losing treatment**", "", "**Winning Treatment**", 
								" ", "Baseline (Control)", " ", "$R^2$", "F", "p", "N")
	return(table)
}

table_Gender <- create_regtable(lm(as.numeric(Gender)-1 ~ Group, df), df)	
table_Age <- create_regtable(lm(Age ~ Group, df), df)
table_EduLevel <- create_regtable(lm(as.numeric(Edu_Level2)+1 ~ Group, df), df)
table_MathConf <- create_regtable(lm(as.numeric(Math_Conf) ~ Group, df), df)
table <- cbind(table_Gender, table_Age, table_EduLevel, table_MathConf)
colnames(table) <- c("Gender", "Age", "Education level", "Math Confidence")
kable(table, align = "r", caption = "Effect of the treatment on covariates")


## @knitr Check_Balanced_Groups_3
hist_Gender <- ggplot(df, aes(Group, fill = Gender))
y0 <- cumsum(prop.table(table(df$Gender))); names(y0) <- NULL
hist_Gender + geom_bar(position="fill") + ylab("Proportion") + 
	scale_fill_discrete(name = "Gender", guide = guide_legend(reverse=TRUE)) + 
	geom_hline(aes(yintercept = y0[1]), colour = "Black", linetype = "dashed", 
				  size = 0.5)




## @knitr Attrition_1
# Some of those 339 records contain missing values (in 2 of the 3 mandatory
	# variables: Score and Score_PRE):
any(is.na(Results$Score)); sum(is.na(Results$Score))
any(is.na(Results$Score_PRE)); sum(is.na(Results$Score_PRE))
any(is.na(Results$Group))

needed_fields <- c("Score_PRE", "Score")

Valid_IDs <- Results$User_ID[apply(Results[, needed_fields], 1, function(x) 
	all(!is.na(x)))]
length(Valid_IDs)
100 - round(100 * length(Valid_IDs) / length(Results$User_ID), 1)
# 307 records (90.6%; attrition rate = 9.4%) has missing values of Score_PRE, 
	# Score, or both.


## @knitr Attrition_2
# We also check for participants that did not follow the instructions and whose
	# results may be inaccurate (and invalid).
# The instructions clearly ordered to use a personal computer, and not a tablet
	# or a smartphone, so those participants who used a device with a small 
	# screen (that limited the speed at which they could enter their answers) are 
	# considered as invalid (attrited).
min(na.omit(Results$Screen_Width))
min(na.omit(Results$Screen_Height))
sort(unique(Results$Screen_Width[Results$User_ID %in% Valid_IDs]))
Results[Results$Screen_Width < 1000 & Results$User_ID %in% Valid_IDs, 
		  c("User_ID", needed_fields, "Screen_Width", "Screen_Height")]
# There were 4 participants that used a 375x667 screen or smaller to finish the
	# experiment. We also consider them as attrited.
Valid_IDs <- Results$User_ID[Results$User_ID %in% Valid_IDs & 
									  	Results$Screen_Width >= 1000]
length(Valid_IDs)
100 - round(100 * length(Valid_IDs) / length(Results$User_ID), 1)
# Now we are left with 303 records (attrition rate = 10.6%).

Results_NM <- Results[Results$User_ID %in% Valid_IDs, ]
row.names(Results_NM) <- NULL
# NM = non-missing
nrow(Results_NM)

df <- Results_NM




## @knitr Graphs
hist_Score <- ggplot(df, aes(Score_PRE))
hist_Score + geom_histogram(binwidth = 5, fill = "White", colour = "Black", 
										  aes(y = ..density..)) + ylab("Frequency") + 
	xlab("1st score") + ggtitle("Histogram of Scores in the 1st test")
hist_Score <- ggplot(df, aes(Score))
hist_Score + geom_histogram(binwidth = 5, fill = "White", colour = "Black", 
									 aes(y = ..density..)) + ylab("Frequency") + 
	xlab("2nd score") + ggtitle("Histogram of Scores in the 2nd test")
hist_Score <- ggplot(df, aes(Score - Score_PRE))
hist_Score + geom_histogram(binwidth = 5, fill = "White", colour = "Black", 
									 aes(y = ..density..)) + ylab("Frequency") + 
	xlab("Score difference") + ggtitle("Histogram of Score differences")
# hist(df$Score_PRE, xlab = "1st score", freq = FALSE, ylim = c(0, 0.05), 
# 	  main = "Histogram of Scores in the 1st test")
# hist(df$Score, xlab = "2nd score", freq = FALSE, ylim = c(0, 0.05), 
# 	  main = "Histogram of Scores in the 2nd test")
# hist(df$Score - df$Score_PRE, xlab = "Score - Score_PRE", freq = FALSE, 
# 	  ylim = c(0, 0.1), main = "Histogram of Score differences")

min_score <- min(c(df$Score, df$Score_PRE))
max_score <- max(c(df$Score, df$Score_PRE))
range_score <- c(floor((min_score - 0.1)/5)*5, ceiling((max_score + 0.1)/5)*5)

density_Score <- ggplot(df, aes(Score_PRE))
density_Score + geom_density(fill = "White") + ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "1st score") + 
	ggtitle("Density plot of Scores in the 1st test") + 
	geom_vline(aes(xintercept = mean(Score_PRE)), colour = "red", 
				  linetype = "dashed", size = 1)
density_Score <- ggplot(df, aes(Score))
density_Score + geom_density(fill = "White") + ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "2nd score") + 
	ggtitle("Density plot of Scores in the 2nd test") + 
	geom_vline(aes(xintercept = mean(Score)), colour = "red", 
				  linetype = "dashed", size = 1)
density_Score <- ggplot(df, aes(Score - Score_PRE))
density_Score + geom_density(fill = "White") + 
	scale_y_continuous(limits = c(0, 0.125), "Frequency") + 
	scale_x_continuous("Score difference") + 
	ggtitle("Density plot of Score differences") + 
	geom_vline(aes(xintercept = mean(Score - Score_PRE)), colour = "red", 
				  linetype = "dashed", size = 1)
# plot(density(df$Score_PRE), xlim = c(0, maximum_score), ylim = c(0, 0.075), 
# 	  yaxp = c(0, 0.075, 15))
# abline(v = mean(df$Score_PRE))
# plot(density(df$Score), xlim = c(0, maximum_score), ylim = c(0, 0.075), 
# 	  yaxp = c(0, 0.075, 15), main = "2nd score")
# abline(v = mean(df$Score))

Mean_Score_PRE_Group <- ddply(df, .(Group), function(.d) 
	data.frame(x = mean(.d$Score_PRE)))
colnames(Mean_Score_PRE_Group)[2] <- "Score_PRE"
Mean_Score_Group <- ddply(df, .(Group), function(.d) 
	data.frame(x = mean(.d$Score)))
colnames(Mean_Score_Group)[2] <- "Score"
Mean_Score_DIFF_Group <- ddply(df, .(Group), function(.d) 
	data.frame(x = mean(.d$Score - .d$Score_PRE)))
colnames(Mean_Score_DIFF_Group)[2] <- "Score_DIFF"
Mean_Score_PRE_Gender <- ddply(df, .(Gender), function(.d) 
	data.frame(x = mean(.d$Score_PRE)))
colnames(Mean_Score_PRE_Gender)[2] <- "Score_PRE"
Mean_Score_Gender <- ddply(df, .(Gender), function(.d) 
	data.frame(x = mean(.d$Score)))
colnames(Mean_Score_Gender)[2] <- "Score"
Mean_Score_DIFF_Gender <- ddply(df, .(Gender), function(.d) 
	data.frame(x = mean(.d$Score - .d$Score_PRE)))
colnames(Mean_Score_DIFF_Gender)[2] <- "Score_DIFF"
Mean_Score_PRE_Group_Gender <- ddply(df, .(Group, Gender), function(.d) 
	data.frame(x = mean(.d$Score_PRE)))
colnames(Mean_Score_PRE_Group_Gender)[3] <- "Score_PRE"
Mean_Score_Group_Gender <- ddply(df, .(Group, Gender), function(.d) 
	data.frame(x = mean(.d$Score)))
colnames(Mean_Score_Group_Gender)[3] <- "Score"
Mean_Score_DIFF_Group_Gender <- ddply(df, .(Group, Gender), function(.d) 
	data.frame(x = mean(.d$Score - .d$Score_PRE)))
colnames(Mean_Score_DIFF_Group_Gender)[3] <- "Score_DIFF"

density_Score <- ggplot(df, aes(Score_PRE))
density_Score + geom_density(fill = "White") + ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "1st score") + 
	ggtitle("Density plot of Scores in the 1st test by Gender") + 
	geom_vline(data = Mean_Score_PRE_Gender, 
				  aes(xintercept = Score_PRE, colour = Gender), 
				  linetype = "dashed", size = 1) + facet_grid(Gender ~ .)
density_Score <- ggplot(df, aes(Score))
density_Score + geom_density(fill = "White") + ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "2nd score") + 
	ggtitle("Density plot of Scores in the 2nd test by Gender") + 
	geom_vline(data = Mean_Score_Gender, 
				  aes(xintercept = Score, colour = Gender), 
				  linetype = "dashed", size = 1) + facet_grid(Gender ~ .)
density_Score <- ggplot(df, aes(Score - Score_PRE))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + scale_x_continuous("Score difference") + 
	ggtitle("Density plot of Score differences by Gender") + 
	geom_vline(data = Mean_Score_DIFF_Gender, 
				  aes(xintercept = Score_DIFF, colour = Gender), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Gender ~ .)

density_Score <- ggplot(df, aes(Score_PRE))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "1st score") + 
	ggtitle("Density plot of Scores in the 1st test by Group") + 
	geom_vline(data = Mean_Score_PRE_Group, 
				  aes(xintercept = Score_PRE, colour = Group), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Group ~ .)
density_Score <- ggplot(df, aes(Score))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "2nd score") + 
	ggtitle("Density plot of Scores in the 2nd test by Group") + 
	geom_vline(data = Mean_Score_Group, 
				  aes(xintercept = Score, colour = Group), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Group ~ .)
density_Score <- ggplot(df, aes(Score - Score_PRE))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + xlab("Score difference") + 
	ggtitle("Density plot of Score differences by Group") + 
	geom_vline(data = Mean_Score_DIFF_Group, 
				  aes(xintercept = Score_DIFF, colour = Group), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Group ~ .)

density_Score <- ggplot(df, aes(Score_PRE))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "1st score") + 
	ggtitle("Density plot of Scores in the 1st test by Group and Gender") + 
	geom_vline(data = Mean_Score_PRE_Group_Gender, 
				  aes(xintercept = Score_PRE, colour = Gender), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Gender ~ Group)
density_Score <- ggplot(df, aes(Score))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + 
	scale_x_continuous(limits = range_score, "2nd score") + 
	ggtitle("Density plot of Scores in the 2nd test by Group and Gender") + 
	geom_vline(data = Mean_Score_Group_Gender, 
				  aes(xintercept = Score, colour = Gender), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Gender ~ Group)
density_Score <- ggplot(df, aes(Score - Score_PRE))
density_Score + geom_density(fill = "White") + 
	ylab("Frequency") + xlab("Score difference") + 
	ggtitle("Density plot of Score differences by Group and Gender") + 
	geom_vline(data = Mean_Score_DIFF_Group_Gender, 
				  aes(xintercept = Score_DIFF, colour = Gender), 
				  linetype = "dashed", size = 1) + 
	facet_grid(Gender ~ Group)

# plot(density(df$Score_PRE[df$Gender == "male"]), xlim = c(0, maximum_score), 
# 	  ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), main = "1st score - Male")
# abline(v = mean(df$Score_PRE[df$Gender == "male"]))
# plot(density(df$Score_PRE[df$Gender == "female"]), xlim = c(0, maximum_score), 
# 	  ylim = c(0, 0.075), yaxp = c(0, 0.075, 15),main = "1st score - Female")
# abline(v = mean(df$Score_PRE[df$Gender == "female"]))
# plot(density(df$Score[df$Gender == "male" & df$Group == "Control"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Control - Male")
# abline(v = mean(df$Score[df$Gender == "male" & df$Group == "Control"]))
# plot(density(df$Score[df$Gender == "female" & df$Group == "Control"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Control - Female")
# abline(v = mean(df$Score[df$Gender == "female" & df$Group == "Control"]))
# plot(density(df$Score[df$Gender == "male" & df$Group == "Loser"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Loser - Male")
# abline(v = mean(df$Score[df$Gender == "male" & df$Group == "Control"]))
# plot(density(df$Score[df$Gender == "female" & df$Group == "Loser"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Loser - Female")
# abline(v = mean(df$Score[df$Gender == "female" & df$Group == "Loser"]))
# plot(density(df$Score[df$Gender == "male" & df$Group == "Winner"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Winner - Male")
# abline(v = mean(df$Score[df$Gender == "male" & df$Group == "Winner"]))
# plot(density(df$Score[df$Gender == "female" & df$Group == "Winner"]), 
# 	  xlim = c(0, maximum_score), ylim = c(0, 0.075), yaxp = c(0, 0.075, 15), 
# 	  main = "2nd score - Winner - Female")
# abline(v = mean(df$Score[df$Gender == "female" & df$Group == "Winner"]))

bar_Score <- ggplot(df, aes(Group, Score_PRE))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge", 
						 fill = "White", colour = "Black") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the 1st mean score by group") + 
	ylab("1st score")
bar_Score <- ggplot(df, aes(Group, Score))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge", 
						 fill = "White", colour = "Black") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the 2nd score mean by group") + 
	ylab("2nd score")
bar_Score <- ggplot(df, aes(Group, Score - Score_PRE))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge", 
						 fill = "White", colour = "Black") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the score difference mean by group") + 
	ylab("Score difference")

bar_Score <- ggplot(df, aes(Group, Score_PRE, fill = Gender))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the 1st score mean by group and gender") + 
	ylab("1st score")
bar_Score <- ggplot(df, aes(Group, Score, fill = Gender))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the 2nd score mean by group and gender") + 
	ylab("2nd score")
bar_Score <- ggplot(df, aes(Group, Score - Score_PRE, fill = Gender))
bar_Score + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + 
	stat_summary(fun.data = mean_cl_normal, geom = "errorbar", 
					 position = position_dodge(width = 0.9), width = 0.2) + 
	ggtitle("Bar chart of the score difference mean by group and gender") + 
	ylab("Score difference")

df_2 <- df[, c("User_ID", "Gender", "Score_PRE", "Score", "Group")]
df_2 <- melt(df_2, id = c("User_ID", "Gender", "Group"), 
				 measured = c("Score_PRE", "Score"))
line_Score <- ggplot(df_2, aes(variable, value, colour = Group))
line_Score + stat_summary(fun.y = mean, geom = "point") + 
	stat_summary(fun.y = mean, geom = "line", aes(group = Group)) + 
	stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + 
	ylab("Score") + scale_x_discrete("", labels = c("1st test", "2nd test")) + 
	ggtitle("Error bar graph of the mean score by group")

line_Score + stat_summary(fun.y = mean, geom = "point") + 
	stat_summary(fun.y = mean, geom = "line", aes(group = Group)) + 
	stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + 
	ylab("Score") + scale_x_discrete("", labels = c("1st test", "2nd test")) + 
	facet_grid(Gender ~ .) + 
	ggtitle("Error bar graph of the mean score by group and gender")

	
	


## @knitr Attrition_3
# We have considered as attrited those participants who left the experiment 
	# (or used a small-screen device)
# What about who did not tried their best in any of the two tests?
	# Those who didn't answer any question at all, or only a few ones
# Similarly, those who answered a large number of questions in one test, and 
	# very few in the other, are likely to have been distracted during the 
	# one with the low score
# To account for this probably non-valid results, we will discard those records
	# (but we will also analyze the whole dataset as well)

boxplot_Score <- ggplot(Results_NM, aes(factor(0), Score_PRE))
boxplot_Score + geom_boxplot() + xlab("") + ylab("1st score") + 
	ggtitle("Boxplot of the 1st score")
boxplot_Score <- ggplot(Results_NM, aes(factor(0), Score))
boxplot_Score + geom_boxplot() + xlab("") + ylab("2nd score") + 
	ggtitle("Boxplot of the 2nd score")
boxplot_Score <- ggplot(Results_NM, aes(factor(0), Score - Score_PRE))
boxplot_Score + geom_boxplot() + xlab("") + ylab("Score difference") + 
	ggtitle("Boxplot of the score difference")

# qtle_PRE <- quantile(Results_NM$Score_PRE, probs = c(0.05, .95))
# qtle <- quantile(Results_NM$Score, probs = c(0.05, .95))
# qtle_DIFF <- quantile(Results_NM$Score - Results_NM$Score_PRE, 
# 							 probs = c(0.05, .95))

# On the one hand, we are going to trim out subjects whose score in any on the 
	# two tests is very low (below the 5th percentile), or whose score 
	# difference is below or above the 5th percentile.
outliers_mdn_Score_PRE <- subset(Results_NM, Score_PRE %in% 
												find_outliers_mdn(Results_NM$Score_PRE, 
																		5, 0), 
											select = c(User_ID, Score_PRE))
outliers_mdn_Score <- subset(Results_NM, Score %in% 
									  	find_outliers_mdn(Results_NM$Score, 5, 0), 
									  select = c(User_ID, Score))
outliers_mdn_Score_DIFF <- subset(Results_NM, (Score - Score_PRE) %in% 
											 	find_outliers_mdn(Results_NM$Score - 
											 								Results_NM$Score_PRE, 5, 5), 
											 select = c(User_ID, Score, Score_PRE))
outliers_mdn <- subset(Results_NM, User_ID %in% outliers_mdn_Score_PRE$User_ID | 
							  	User_ID %in% outliers_mdn_Score$User_ID | 
							  	User_ID %in% outliers_mdn_Score_DIFF$User_ID, 
							  select = c(User_ID, Score, Score_PRE))
outliers_mdn_Score_DIFF$Score_DIFF <- outliers_mdn_Score_DIFF$Score - 
	outliers_mdn_Score_DIFF$Score_PRE
outliers_mdn$Score_DIFF <- outliers_mdn$Score - outliers_mdn$Score_PRE

IDs_trm_mdn <- Results_NM$User_ID[!(Results_NM$User_ID %in% 
													outliers_mdn$User_ID)]

length(sort(Results_NM$Score[Results_NM$User_ID %in% IDs_trm_mdn]))
# This way we discard 43 observations.

# On the other hand, we are going to assume normality and trim out subjects 
# whose score in any on the two tests is very low (probability < 2.5%) or
# whose absolute score difference is very big (probability < 5% in both 
# directions).
outliers_mn_Score_PRE <- subset(Results_NM, Score_PRE %in% 
										  	find_outliers_mn(Results_NM$Score_PRE, 2.5, 0), 
										  select = c(User_ID, Score_PRE))
outliers_mn_Score <- subset(Results_NM, Score %in% 
									 	find_outliers_mn(Results_NM$Score, 2.5, 0), 
									 select = c(User_ID, Score))
outliers_mn_Score_DIFF <- subset(Results_NM, (Score - Score_PRE) %in% 
												find_outliers_mn(Results_NM$Score - 
																	  	Results_NM$Score_PRE, 5, 
																	  5), 
											select = c(User_ID, Score, Score_PRE))
outliers_mn <- subset(Results_NM, User_ID %in% outliers_mn_Score_PRE$User_ID | 
							 	User_ID %in% outliers_mn_Score$User_ID | 
							 	User_ID %in% outliers_mn_Score_DIFF$User_ID, 
							 select = c(User_ID, Score, Score_PRE))
outliers_mn_Score_DIFF$Score_DIFF <- outliers_mn_Score_DIFF$Score - 
	outliers_mn_Score_DIFF$Score_PRE
outliers_mn$Score_DIFF <- outliers_mn$Score - outliers_mn$Score_PRE

IDs_trm_mn <- Results_NM$User_ID[!(Results_NM$User_ID %in% 
											  	outliers_mn$User_ID)]

length(sort(Results_NM$Score[Results_NM$User_ID %in% IDs_trm_mn]))
# This way we discard 41 observations.

outliers_mn[which(!(outliers_mn$User_ID %in% outliers_mdn$User_ID)), ]
outliers_mdn[which(!(outliers_mdn$User_ID %in% outliers_mn$User_ID)), ]
# Both methods discard almost the same observations, so we will just use one of 
# them (the one which discards abnormal values with respect to the means).

Results_NM_mdn <- Results_NM[-which(Results_NM$User_ID %in% 
													outliers_mdn$User_ID), ]
row.names(Results_NM_mdn) <- NULL

Results_NM_mn <- Results_NM[-which(Results_NM$User_ID %in% 
											  	outliers_mn$User_ID), ]
row.names(Results_NM_mn) <- NULL

# Now the distributions are closer to normality
by(Results_NM$Score_PRE, Results_NM$Group, shapiro.test)
by(Results_NM_mn$Score_PRE, Results_NM_mn$Group, shapiro.test)

by(Results_NM$Score, Results_NM$Group, shapiro.test)
by(Results_NM_mn$Score, Results_NM_mn$Group, shapiro.test)

by(Results_NM$Score - Results_NM$Score_PRE, Results_NM$Group, shapiro.test)
by(Results_NM_mn$Score - Results_NM_mn$Score_PRE, Results_NM_mn$Group, 
	shapiro.test)




## @knitr Descriptive_Statistics
descriptive <- function(df) {
	winner <- paste0(formatC(100*mean(df$Group == "Winner"), digits = 1, 
									 format = "f"), "%")
	loser <- paste0(formatC(100*mean(df$Group == "Loser"), digits = 1, 
									format = "f"), "%")
	score1 <- formatC(mean(df$Score_PRE), digits = 1, format = "f")
	score2 <- formatC(mean(df$Score), digits = 1, format = "f")
	male <- paste0(formatC(100*mean(df$Gender == "male"), digits = 1, 
								  format = "f"), "%")
	age <- formatC(mean(df$Age), digits = 1, format = "f")
	Edu_Level2 <- factor(df$Edu_Level, levels = c(1:9))
	edulevel <- formatC(mean(as.numeric(Edu_Level2)), digits = 1, format = "f")
	N <- nrow(df)
	table <- matrix(c(winner, loser, score1, score2, male, age, edulevel, N), 
						 ncol = 1)
	rownames(table) <- c("Assigned to winning condition", 
								"Assigned to losing condition", 
								"Score in the 1st test (mean)", 
								"Score in the 2nd test (mean)", "Male", "Age (mean)", 
								"Education level (mean)", "N (sample size)")
	return(table)
}

# table <- mapply(cbind, lapply(list(Results_NM, Results_NM_mdn), 
# 										descriptive))
table <- cbind(descriptive(Results_NM), descriptive(Results_NM_mdn))
colnames(table) <- c("w/o missing values", "w/o abnormal values")
kable(table, align = "r", 
		caption = "Descriptive statistics of the sample under study")




## @knitr Attrition_Analysis_1
df <- Results_NM

## @knitr Attrition_Analysis_2
df_whole <- Results

NonObserved_IDs <- df_whole$User_ID[(!df_whole$User_ID %in% df$User_ID)]
df_whole$observed <- ifelse(df_whole$User_ID %in% NonObserved_IDs, 0, 1) 

df_whole$probobs <- glm(observed ~ (Group*Age) + (Group*Gender) + 
									(Group*Edu_Level) + (Group*Math_Conf), 
								data = df_whole, 
								family = binomial(link = "logit"))$fitted
# Compare distributions of predicted probabilities across experimental 
	# conditions.
# Make sure that there are no zero predicted probabilities in either condition.
by(df_whole$probobs, df_whole$Group, summary)

# Generate weights: inverse of predicted probability of being observed.
df_whole$wt <- 1 / df_whole$probobs

# Restrict analysis to observed subjects.
df_whole$sel_valid <- df_whole$observed == 1
table(df_whole$sel_valid)

# Coefficients for unweighted regression (restricting analysis to observed 
	# subjects).
lm(Score ~ Group, subset = sel_valid, df_whole)$coefficients
cl(lm(Score ~ Group, subset = sel_valid, df_whole), 
	df_whole$Cluster_ID[df_whole$sel_valid == 1])
# Coefficients for IPW regression (restricting analysis to observed subjects).
lm(Score ~ Group, weights = wt, subset=sel_valid, df_whole)$coefficients
cl(lm(Score ~ Group, weights = wt, subset=sel_valid, df_whole), 
	df_whole$Cluster_ID[df_whole$sel_valid == 1])

cl(lm(Score ~ Group*Gender, subset = sel_valid, df_whole), 
	df_whole$Cluster_ID[df_whole$sel_valid == 1])
cl(lm(Score ~ Group*Gender, weights = wt, subset=sel_valid, df_whole), 
	df_whole$Cluster_ID[df_whole$sel_valid == 1])

# Regression estimates predicting missingness as a function of covariates, by 
	# experimental Group:
df_whole$missing <- 1 - df_whole$observed

attr_rates <- paste0(formatC(as.numeric(100*by(df_whole$missing, df_whole$Group, 
															 mean)), digits = 1, 
									 format = "f", drop0trailing = FALSE), "%")
ANOVA <- summary.lm(aov(missing ~ Group, df_whole))

#Regression for Control Group
summary(lm(missing ~ Gender + Edu_Level + Math_Conf, 
			  subset = (Group == "Control"), df_whole))$coefficients
cl(lm(missing ~ Gender + Edu_Level + Math_Conf, 
			  subset = (Group == "Control"), df_whole), 
	df_whole$Cluster_ID[df_whole$Group == "Control"])
#Regression for Loser Group
summary(lm(missing ~ Gender + Edu_Level + Math_Conf, 
			  subset = (Group == "Loser"), df_whole))$coefficients
cl(lm(missing ~ Gender + Edu_Level + Math_Conf, 
		subset = (Group == "Loser"), df_whole), 
	df_whole$Cluster_ID[df_whole$Group == "Loser"])
#Regression for Treatment Group
summary(lm(missing ~ Gender + Edu_Level + Math_Conf, 
			  subset = (Group == "Winner"), df_whole))$coefficients
cl(lm(missing ~ Age + Gender + Edu_Level + Math_Conf, 
		subset = (Group == "Winner"), df_whole), 
	df_whole$Cluster_ID[df_whole$Group == "Winner"])

Groups <- c("Control", "Loser", "Winner")
table1 <- Map(function(x) 
	cl(lm(missing ~ Age + Gender + Edu_Level + Math_Conf, subset = Group == x, 
			df_whole), df_whole$Cluster_ID[df_whole$Group == x]), Groups)
# table1 <- Map(function(x) 
# 	summary(lm(missing ~ Age + Gender + Edu_Level + Math_Conf, subset =  Group == x, 
# 				  df_whole))$coefficients, Groups)
table2 <- lapply(table1, function(x) 
	apply(x, 1, function(y) paste0(formatC(y[1], digits = 3, format = "f", 
														drop0trailing = FALSE), 
											 sig_stars(y[4]), " (", 
											 formatC(y[2], digits = 3, 
											 		  format = "f", drop0trailing = FALSE), 
											 ")")))
cats <- c(row.names(table1[[1]])[1:3], 
			 unlist(lapply(levels(df_whole$Edu_Level)[
			 	2:length(levels(df_whole$Edu_Level))], function(x) 
			 		paste0("Edu_Level", x))), 
			 unlist(lapply(levels(df_whole$Math_Conf)[
			 	2:length(levels(df_whole$Math_Conf))], function(x) 
			 		paste0("Math_Conf", x))))
for (i in 1:3) {
	for (j in cats[!(cats %in% names(table2[[Groups[i]]]))]) {
		table2[[i]][j] <- ""
	}
	table2[[i]] <- table2[[i]][cats]
}
table3 <- matrix(unlist(table2), byrow = FALSE, nrow = 14, ncol = 3)
rownames(table3) <- c("Baseline", "Age", "Female", "Education level 7", 
							 "Education level 8", "Education level 9", 
							 "Education level 5", "Education level 4", 
							 "Education level 3", "Education level 2", 
							 "Math confidence 2", "Math confidence 3", 
							 "Math confidence 4", "Math confidence 5")
colnames(table3) <- Groups


## @knitr Attrition_Analysis_3
kable(table3, align = "r")


## @knitr Attrition_Analysis_4
create_regtable3 <- function(model, df) {
	model2 <- summary(model)
	model3 <- cl(model, df$Cluster_ID[df$sel_valid == 1])
	estimate <- unlist(lapply(c(1:6), function(x) 
		paste0(formatC(model3[x, 1], digits = 3, format = "f", 
							drop0trailing = FALSE), sig_stars(model3[x, 4]))))
	SE <- unlist(lapply(c(1:6), function(x) 
		paste0("(", formatC(model3[x, 2], digits = 3, format = "f", 
								  drop0trailing = FALSE), ")  ")))
	N <- paste0(nrow(df), "   ")
	attrited <- paste0(formatC(100*mean(df$sel_valid==0), digits = 1, 
										format = "f", drop0trailing = FALSE) , "%   ")
	R2 <- paste0(formatC(model2$r.squared, digits = 3, format = "f", 
								drop0trailing = FALSE), "   ")
	Fstatistic <- paste0(formatC(model2$fstatistic[1], digits = 3, format = "f", 
										  drop0trailing = FALSE), "   ")
	pvalue <- paste0(formatC(1 - pf(model2$fstatistic[1], 2, 300), digits = 3, 
									 format = "f", drop0trailing = FALSE), "   ")
	values <- c(t(matrix(c(estimate, SE), ncol = 2)), R2, Fstatistic, pvalue, 
					N, attrited)
	table <- data.frame(variables = c(c(1:6), "Covariate", "Covariate (SE)", 
												 "L:Covariate", "L:Covariate (SE)", 
												 "W:Covariate", "W:Covariate (SE)", 
												 c(13:17)), values)
	return(table)
}

reg1 <- create_regtable3(lm(Score ~ Group * Score_PRE, weights = wt, 
									 subset=sel_valid, df_whole), df_whole)
reg1[, 1] <- gsub("Covariate", "Score_PRE", reg1[, 1])
reg2 <- create_regtable3(lm(Score ~ Group * Gender, weights = wt, 
									 subset=sel_valid, df_whole), df_whole)
reg2[, 1] <- gsub("Covariate", "Gender", reg2[, 1])
reg3 <- create_regtable3(lm(Score ~ Group * Age, weights = wt, 
									 subset=sel_valid, df_whole), df_whole)
reg3[, 1] <- gsub("Covariate", "Age", reg3[, 1])
table <- merge(merge(reg1, reg2, by = "variables", all = TRUE), reg3, 
					by = "variables", all = TRUE)
rownames(table) <- table[, 1]
table <- table[c(1, 2:6, "Score_PRE", "Score_PRE (SE)", "L:Score_PRE", 
					  "L:Score_PRE (SE)", "W:Score_PRE", "W:Score_PRE (SE)", 
					  "Gender", "Gender (SE)", "L:Gender", "L:Gender (SE)", 
					  "W:Gender", "W:Gender (SE)", "Age", "Age (SE)", 
					  "L:Age", "L:Age (SE)", "W:Age", "W:Age (SE)", 13:17), ]
table$variables <- c("**Baseline**", "", "**Losing Treatment**", "", 
							"**Winning Treatment**", "", 
							"**(1 point more in) 1st score**", "", 
							"**Losing:1st score**", "", "**Winning:1st score**", "", 
							"**Female**", "", "**Losing:Female**", "", 
							"**Winning:Female**", "", 
							"**Age (per year)**", "", "**Losing:Age**", "", 
							"**Winning:Age**", "", "$R^2$", "F", "p", "N", 
							"Attrition Rate")
rownames(table) <- NULL
table <- sapply(table, as.character)
table[is.na(table)] <- ""
colnames(table) <- c("", "1st score", "Gender", "Age")
kable(table, align = "r")

# a <- by(df_whole$observed, df_whole$Group, mean)
# b <- by(df$Score, df$Group, mean)
# lapply(c(2:3), function(x) 
# 	(a[x]*b[x] + (1 - a[x])*60) - (a[1]*b[1] + (1 - a[1])*0))
# lapply(c(2:3), function(x) 
# 	(a[x]*b[x] + (1 - a[x])*0) - (a[1]*b[1] + (1 - a[1])*60))
# c <- cl(lm(Score ~ Group, df), df$Cluster_ID)
# lapply(c(2:3), function(x) 
# 	c[x, 1] - qnorm(1-.05/2/2)*c[x, 2]) 
# lapply(c(2:3), function(x) 
# 	c[x, 1] + qnorm(1-.05/2/2)*c[x, 2]) 



## @knitr Attrition_Analysis_5
df <- Results_NM_mdn

# cl(lm(missing ~ Group, df_whole), df_whole$Cluster_ID)
# cl(lm(missing ~ Group, df_whole[df_whole$Group != "Loser", ]), 
# 	df_whole$Cluster_ID[df_whole$Group != "Loser"])




## @knitr Pre-Analysis_1
same_cluster <- Results_NM$Cluster_ID[Results_NM$Cluster_ID %in% 
					 	Results_NM$Cluster_ID[duplicated(Results_NM$Cluster_ID)]]
#Some comprobations prior to the final analysis, using the dataset where
# abnormal values (compared to the means) have been discarded:
model_noCovariates <- lm(Score ~ Group, df)
# Robust Standard Errors
# RSEs(model_noCovariates)

# Cluster Standard Errors
# cl(model_noCovariates, df$Cluster_ID)

# Since we performed clustered allocation of subjects, we have to perform the 
# corresponding correction of the SEs.
# But Clustering was made by IP address, so most of the clusters contain just 
# one individual, and hence CSEs are very similar to RSEs
# (The assumption of homoskedasticity is met---variance across Groups is very 
# similar, so RSEs are also quite similar to the SEs directly reporte by 
# the 'lm' function.)
# summary(model_noCovariates)




## @knitr Pre-Analysis_2
mean_score <- rbind(calculate_mean(Results_NM, "Score_PRE"), 
						  calculate_mean(Results_NM, "Score"), 
						  calculate_mean(Results_NM, "Score", "Score_PRE"))
rownames(mean_score) <- c("Mean score on 1st test", "Mean score on 2nd test", 
								  "Mean score difference")
kable(mean_score, align = "r", caption = "Mean score by group")


## @knitr Pre-Analysis_3
df <- Results_NM
plot1 <- plot_density(df)
plot2 <- plot_density_DIFF(df)
grid.arrange(plot1, plot2, ncol=2)


## @knitr Pre-Analysis_4
plot_line(df)

## @knitr Pre-Analysis_5
mean_score <- rbind(calculate_mean(Results_NM_mdn, "Score_PRE"), 
						  calculate_mean(Results_NM_mdn, "Score"), 
						  calculate_mean(Results_NM_mdn, "Score", "Score_PRE"))
rownames(mean_score) <- c("Mean score on 1st test", "Mean score on 2nd test", 
								  "Mean score difference")
kable(mean_score, align = "r", 
		caption = "Mean score by group when abnormal scores are trimmed out")


## @knitr Pre-Analysis_6
df <- Results_NM_mdn
plot1 <- plot_density(df)
plot2 <- plot_density_DIFF(df)
grid.arrange(plot1, plot2, ncol=2)


## @knitr Pre-Analysis_7
plot_line(df)




## @knitr Analysis_1
reg <- create_regtable(lm(Score ~ Group, Results_NM), Results_NM)
reg_mdn <- create_regtable(lm(Score ~ Group, Results_NM_mdn), 
									 Results_NM_mdn)
table <- cbind(reg, reg_mdn)
colnames(table) <- c("w/o missing values", "w/o abnormal values")
kable(table, align = "r", 
		caption = "Effect of the treatment on the score")


## @knitr Analysis_2
reg <- create_regtable(lm(Score - Score_PRE ~ Group, Results_NM), Results_NM)
reg_mdn <- create_regtable(lm(Score - Score_PRE ~ Group, Results_NM_mdn), 
									Results_NM_mdn)
table <- cbind(reg, reg_mdn)
colnames(table) <- c("w/o missing values", "w/o abnormal values")
kable(table, align = "r", 
		caption = "Effect of the treatment on the score difference")


## @knitr Analysis_Hetero_1
create_regtable2 <- function(model, df) {
	model2 <- summary(model)
	model3 <- cl(model, df$Cluster_ID)
	estimate <- unlist(lapply(c(1:6), function(x) 
		paste0(formatC(model3[x, 1], digits = 3, format = "f", 
							drop0trailing = FALSE), sig_stars(model3[x, 4]))))
	SE <- unlist(lapply(c(1:6), function(x) 
		paste0("(", formatC(model3[x, 2], digits = 3, format = "f", 
								  drop0trailing = FALSE), ")  ")))
	N <- paste0(nrow(df), "   ")
	R2 <- paste0(formatC(model2$r.squared, digits = 3, format = "f", 
								drop0trailing = FALSE), "   ")
	Fstatistic <- paste0(formatC(model2$fstatistic[1], digits = 3, format = "f", 
										  drop0trailing = FALSE), "   ")
	pvalue <- paste0(formatC(1 - pf(model2$fstatistic[1], 2, 300), digits = 3, 
									 format = "f", drop0trailing = FALSE), "   ")
	values <- c(t(matrix(c(estimate, SE), ncol = 2)), R2, Fstatistic, pvalue, N)
	table <- data.frame(variables = c(c(1:6), "Covariate", "Covariate (SE)", 
												 "L:Covariate", "L:Covariate (SE)", 
												 "W:Covariate", "W:Covariate (SE)", 
												 c(13:16)), values)
	return(table)
}

reg1 <- create_regtable2(lm(Score ~ Group * Score_PRE, Results_NM), Results_NM)
reg1[, 1] <- gsub("Covariate", "Score_PRE", reg1[, 1])
reg2 <- create_regtable2(lm(Score ~ Group * Gender, Results_NM), Results_NM)
reg2[, 1] <- gsub("Covariate", "Gender", reg2[, 1])
reg3 <- create_regtable2(lm(Score ~ Group * Age, Results_NM), Results_NM)
reg3[, 1] <- gsub("Covariate", "Age", reg3[, 1])
table <- merge(merge(reg1, reg2, by = "variables", all = TRUE), reg3, 
					by = "variables", all = TRUE)
rownames(table) <- table[, 1]
table <- table[c(1, 2:6, "Score_PRE", "Score_PRE (SE)", "L:Score_PRE", 
					  "L:Score_PRE (SE)", "W:Score_PRE", "W:Score_PRE (SE)", 
					  "Gender", "Gender (SE)", "L:Gender", "L:Gender (SE)", 
					  "W:Gender", "W:Gender (SE)", "Age", "Age (SE)", 
					  "L:Age", "L:Age (SE)", "W:Age", "W:Age (SE)", 13:16), ]
table$variables <- c("**Baseline**", "", "**Losing Treatment**", "", 
							"**Winning Treatment**", "", 
							"**(1 point more in) 1st score**", "", 
							"**Losing:1st score**", "", "**Winning:1st score**", "", 
							"**Female**", "", "**Losing:Female**", "", 
							"**Winning:Female**", "", 
							"**Age (per year)**", "", "**Losing:Age**", "", 
							"**Winning:Age**", "", "$R^2$", "F", "p", "N")
rownames(table) <- NULL
table <- sapply(table, as.character)
table[is.na(table)] <- ""
colnames(table) <- c("", "1st score", "Gender", "Age")
kable(table, align = "r", 
		caption = "Heterogeneous effects")


## @knitr Analysis_Hetero_2
scatter <- ggplot(Results_NM, aes(Age, Score, colour = Group))
scatter + geom_point() + geom_smooth(method = "lm", aes(fill = Group), 
												 alpha = 0.2)


## @knitr Analysis_Hetero_3
reg1 <- create_regtable2(lm(Score ~ Group * Score_PRE, Results_NM_mdn), 
								 Results_NM_mdn)
reg1[, 1] <- gsub("Covariate", "Score_PRE", reg1[, 1])
reg2 <- create_regtable2(lm(Score ~ Group * Gender, Results_NM_mdn), 
								 Results_NM_mdn)
reg2[, 1] <- gsub("Covariate", "Gender", reg2[, 1])
reg3 <- create_regtable2(lm(Score ~ Group * Age, Results_NM_mdn), 
								 Results_NM_mdn)
reg3[, 1] <- gsub("Covariate", "Age", reg3[, 1])
table <- merge(merge(reg1, reg2, by = "variables", all = TRUE), reg3, 
					by = "variables", all = TRUE)
rownames(table) <- table[, 1]
table <- table[c(1, 2:6, "Score_PRE", "Score_PRE (SE)", "L:Score_PRE", 
					  "L:Score_PRE (SE)", "W:Score_PRE", "W:Score_PRE (SE)", 
					  "Gender", "Gender (SE)", "L:Gender", "L:Gender (SE)", 
					  "W:Gender", "W:Gender (SE)", "Age", "Age (SE)", 
					  "L:Age", "L:Age (SE)", "W:Age", "W:Age (SE)", 13:16), ]
table$variables <- c("**Baseline**", "", "**Losing Treatment**", "", 
							"**Winning Treatment**", "", 
							"**(1 point more in) 1st score**", "", 
							"**Losing:1st score**", "", "**Winning:1st score**", "", 
							"**Female**", "", "**Losing:Female**", "", 
							"**Winning:Female**", "", 
							"**Age (per year)**", "", "**Losing:Age**", "", 
							"**Winning:Age**", "", "$R^2$", "F", "p", "N")
rownames(table) <- NULL
table <- sapply(table, as.character)
table[is.na(table)] <- ""
colnames(table) <- c("", "1st score", "Gender", "Age")
kable(table, align = "r", 
		caption = "Heterogeneous effects when abnormal scores are trimmed out")





## @knitr test
# the difference-in-differences.
# We prefer this method because the 1st score does not perfectly predict the
	# 2nd score (even in the absence of treatment).
# Anyway, the estimates of the treatments' effects are somewhat similar using 
	# both methods, and both have the effect of reducing the uncertainty of our
	# estimates.
cl(model_noCovariates, df$Cluster_ID)
model_ScorePRE <- lm(Score ~ Score_PRE + Group, df)
cl(model_ScorePRE, df$Cluster_ID)
model_ScoreDIFF <- lm(Score - Score_PRE ~ Group, df)
cl(model_ScoreDIFF, df$Cluster_ID)
# We will use the score in the 1st test as a covariate, rather than analyzing
	# the difference-in-differences.
# We prefer this method because the 1st score does not perfectly predict the
	# 2nd score (even in the absence of treatment).
# Anyway, the estimates of the treatments' effects are somewhat similar using 
	# both methods, and both have the effect of reducing the uncertainty of our
	# estimates.
cl(model_noCovariates, df$Cluster_ID)
model_ScorePRE <- lm(Score ~ Score_PRE + Group, df)
cl(model_ScorePRE, df$Cluster_ID)
model_ScoreDIFF <- lm(Score - Score_PRE ~ Group, df)
cl(model_ScoreDIFF, df$Cluster_ID)
 rather than analyzing
	# the difference-in-differences.
# We prefer this method because the 1st score does not perfectly predict the
	# 2nd score (even in the absence of treatment).
# Anyway, the estimates of the treatments' effects are somewhat similar using 
	# both methods, and both have the effect of reducing the uncertainty of our
	# estimates.
cl(model_noCovariates, df$Cluster_ID)
model_ScorePRE <- lm(Score ~ Score_PRE + Group, df)
cl(model_ScorePRE, df$Cluster_ID)
model_ScoreDIFF <- lm(Score - Score_PRE ~ Group, df)
cl(model_ScoreDIFF, df$Cluster_ID)
# boot.ci(model_ScorePRE_Gender)
c <- bootcoefs(lmrob(Score ~ Group * Age, df))
plot(c)
# This robust version of the regression shows that the positive effect of being 
	# compared to someone who's performing better is near statistically 
	# significant among Women

c <- bootcoefs(lmrob(Score ~ Group, Results_NM))
plot(c)
Results_NM$Score_DIFF <- Results_NM$Score - Results_NM$Score_PRE
c <- bootcoefs(lmrob(Score_DIFF ~ Group, Results_NM))
plot(c)




## @knitr Statistical_Power_1
df <- Results_NM

x<-with(df, by(df, Group, function(x) 
	var(x$Score - x$Score_PRE)))
y<-with(df, by(df, Group, function(x) 
	var(x$Score)))
# y/x
a <- with(df, by(df, Group, function(x) 
	sd(x$Score - x$Score_PRE) / mean(x$Score - x$Score_PRE)))
b <- with(df, by(df, Group, function(x) sd(x$Score) / mean(x$Score)))
# a/b
c <- with(df, by(df, Group, function(x) 
	sd(x$Score - x$Score_PRE) / mean(x$Score - x$Score_PRE)))
d <- with(df, by(df, Group, function(x) 
	sd(x$Score) / mean(x$Score)))
# c/d

d_Loser <- cohen.d(df$Score[df$Group == "Loser"], 
						 df$Score[df$Group == "Control"])
d_Winner <- cohen.d(df$Score[df$Group == "Winner"], 
						 df$Score[df$Group == "Control"])

pwr.t.test(n = nrow(df)/3, d = d_Loser$estimate, sig.level = 0.05/2, 
			  power = NULL)
pwr.t.test(n = NULL, d = d_Loser$estimate, sig.level = 0.05/2, power = 0.8)
pwr.t.test(n = nrow(df)/3, d = d_Winner$estimate, sig.level = 0.05/2, 
			  power = NULL)
pwr.t.test(n = NULL, d = d_Winner$estimate, sig.level = 0.05/2, 
			  power = 0.8)

(model <- summary(lm(Score ~ Group, df)))
pwr <- pwr.f2.test(u = 1, v = nrow(df)-2, f2 = model$r.squared / 
					(1-model$r.squared), sig.level = 0.05/2)
N <- pwr.f2.test(u = 1, v = NULL, f2 = model$r.squared / 
					(1-model$r.squared), sig.level = 0.05/2, power = 0.8)

# estimate_power(df, nrow(df))
# estimate_power(df, 3e3)

# t <- data.frame(Sample_Size = seq(floor(nrow(df)/300)*300, 3e3, 300), 
# 					 Power = unlist(lapply(seq(floor(nrow(df)/300)*300, 3e3, 300), 
# 					 							 function(x) estimate_power(df, x))))
# kable(t)




## @knitr Statistical_Power_2
df <- Results_NM_mdn

x<-with(df, by(df, Group, function(x) 
	var(x$Score - x$Score_PRE)))
y<-with(df, by(df, Group, function(x) 
	var(x$Score)))
# y/x
a <- with(df, by(df, Group, function(x) 
	sd(x$Score - x$Score_PRE) / mean(x$Score - x$Score_PRE)))
b <- with(df, by(df, Group, function(x) sd(x$Score) / mean(x$Score)))
# a/b
c <- with(df, by(df, Group, function(x) 
	sd(x$Score - x$Score_PRE) / mean(x$Score - x$Score_PRE)))
d <- with(df, by(df, Group, function(x) 
	sd(x$Score) / mean(x$Score)))
# c/d

d_Loser <- cohen.d(df$Score[df$Group == "Loser"], 
						 df$Score[df$Group == "Control"])
d_Winner <- cohen.d(df$Score[df$Group == "Winner"], 
						  df$Score[df$Group == "Control"])

pwr.t.test(n = nrow(df)/3, d = d_Loser$estimate, sig.level = 0.05/2, 
			  power = NULL)
pwr.t.test(n = NULL, d = d_Loser$estimate, sig.level = 0.05/2, power = 0.8)
pwr.t.test(n = nrow(df)/3, d = d_Winner$estimate, sig.level = 0.05/2, 
			  power = NULL)
pwr.t.test(n = NULL, d = d_Winner$estimate, sig.level = 0.05/2, 
			  power = 0.8)

(model <- summary(lm(Score ~ Group, df)))
pwr <- pwr.f2.test(u = 1, v = nrow(df)-2, f2 = model$r.squared / 
								(1-model$r.squared), sig.level = 0.05/2)
N <- pwr.f2.test(u = 1, v = NULL, f2 = model$r.squared / 
					(1-model$r.squared), sig.level = 0.05/2, power = 0.8)

# estimate_power(df, nrow(df))
# estimate_power(df, 7.4e3)

# t <- data.frame(Sample_Size = seq(floor(nrow(df)/300)*300, 3e3, 300), 
# 					 Power = unlist(lapply(seq(floor(nrow(df)/300)*300, 3e3, 300), 
# 					 							 function(x) estimate_power(df, x))))
# kable(t)
