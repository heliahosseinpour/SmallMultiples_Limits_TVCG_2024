source(file = "Functions.R")
library(tidyverse)


df <- read.csv(file = "raw_data/Lookup1 Norming_October 16, 2023_14.50.csv")
#Identify1 - remove unwanted columns

OrderKey <- df[-1:-2, -c(1:8, 10:27) ] %>% 
  dplyr::mutate(task = "Identify1") %>% #Task Id column
  dplyr::select(ResponseId, FL_60_DO) %>% 
  mutate(numbers = stringr::str_extract_all(FL_60_DO, "\\d+")) %>% 
  tidyr::separate(col = numbers,into = as.character(1:18),sep = ",") %>% 
  dplyr::mutate(across(`1`:`18`, function(x) str_extract(x, "\\d+"))) %>% 
  tidyr::pivot_longer(cols = `1`:`18`,
                      names_to = "Order", #Where the column names go
                      values_to = "Frame_Quantity", #Where the values in thos columns go
                      values_drop_na = TRUE) %>% 
  dplyr::select(ResponseId,Frame_Quantity = as.character("Frame_Quantity"),Order)

Identify1_LongV3 <- read.csv(file = "raw_data/Experiment1_Identify1_Raw.csv")[ -c(1:8, 10:26) ] %>% 
  dplyr::mutate(task = "Identify1") %>% #Task Id column
  tidyr::pivot_longer( #Pivot to long format for ID1 Trials and Frames
    cols = Look1_2GraphsFill1_1:Look1_70GraphsFill20_70, #Columns of all the Frames
    names_to = "trial_info", #Where the column names go
    values_to = "Response", #Where the values in thos columns go
    values_drop_na = TRUE)%>% #Drop NAs
  mutate(ResponseCode = as.numeric(ifelse(Response   == "On", "1", "0"))) %>% #Make On and Off numeric (0,1)
  tidyr::separate(trial_info, into = c("task1","trialinfo2","Individual_FramesNum"), #names column split into 3
                  sep = "_",remove = FALSE) %>%  #split column to make Individual_FramesNum
  mutate(Frame_Quantity = str_remove(trialinfo2, "Graphs.*"), #Create FrameQuantity
         Trial = str_remove(trialinfo2,'.*Fill'), #Create Trial
         Trial = str_remove(Trial, '\\.')) %>% #Remove extra periods from the end of trial
  select(-trialinfo2,-task1, -Gender_4_TEXT,-Educate_9_TEXT, #Remove extra columns
         -trial_info, -Race_8_TEXT,-Response) %>% 
  relocate(c(Individual_FramesNum,ResponseCode) , .after = Trial) %>% #Put IFN and RC after Trial
  mutate(across(.cols = everything(),.fns = as.character)) #Make everything a character for joining purposes


Identify1_LongV3 <- left_join(Identify1_LongV3,OrderKey, by = c("ResponseId", "Frame_Quantity"))


load(file = "stimuli/fullkey.RData") #Load the LIST OF LISTS (aka, the masterkey)
#Fix company name to uniform name scheme (in this case upper case letters)
listOlist_dats[[1]][[1]] <- data.frame(Date = listOlist_dats[[1]][[1]]$Date,
                                       my_vector = listOlist_dats[[1]][[1]]$my_vector,
                                       company_name = toupper(c(rep(letters[1],13),rep(letters[2],13),rep(letters[3],13),rep(letters[4],13),
                                                                rep(letters[5],13),rep(letters[6],13))))

#Get max values for each frame
frameMaxes <- lapply(listOlist_dats, function(x) lapply(x, function(y) tapply(y[,2],INDEX = y[,3],FUN = max)))
#Subset 2 frames from A,B,C,D,E,F to just A,F
frameMaxes[[1]] <-  lapply(frameMaxes[[1]], function(x) x[c("A", "F")])

#Rename trials from letters (eg A,F) to numbers
for (i in seq_along(frameMaxes)) { #Iterate through the list of 18 (1st level - frame#) in FrameMaxes
  frameMaxes[[i]] <- lapply(frameMaxes[[i]], function(x) { #Iterate through L20 (2nd level - trials)
    names(x) <- seq(1, length(x), 1) #Name each frame depending on frame#
    return(x)
  })
}

# Might come in handy later...
Maxes <- lapply(frameMaxes, function(x) lapply(x, max)) #Get the Max (ie winner) value for each trial

MaxPos <- lapply(frameMaxes, function(x) { #Get the Max (ie winner) FRAME for each trial
  lapply(x, function(y) {
    max_index <- which.max(y)
    names(y)[max_index]
  })
})

names(MaxPos) <- unique(Identify1_LongV3$Frame_Quantity) #Match the Frame Quantities to the answer key

C_ansr <- Identify1_LongV3 %>%  # Get only the frames that were selected by a participant.
  filter(ResponseCode == 1) 

# create an empty data.frame to store the results
key_Identify_df <- data.frame(Frame_Quantity = integer(),
                              Trial = integer(),
                              Correct_FramesNum = integer(),
                              stringsAsFactors = FALSE)

# loop through each list in the main list
for (i in seq_along(MaxPos)) {
  # loop through each element in the list
  for (j in seq_along(MaxPos[[i]])) {
    # extract the values and add them to the data.frame
    key_Identify_df <- rbind(key_Identify_df, data.frame(Frame_Quantity = names(MaxPos)[i],
                                                         Trial = as.character(j),
                                                         Correct_FramesNum = MaxPos[[i]][[j]],
                                                         stringsAsFactors = FALSE))
  }
}

dfjoined <- left_join(C_ansr, key_Identify_df, by = c("Frame_Quantity","Trial"))

Exp1_ID1_df <- dfjoined


#### Order Analysis


dfaccID1 <- as_tibble(Exp1_ID1_df %>% 
                        # dplyr::filter(ResponseId %ni% c("R_dc1qAAvQQo7m133",
                        #               "R_2QVBIR0ja3iFuUL",
                        #               "R_3R8I8VWOK8N3ZIf",
                        #               "R_cNLo7w2zn57UVhL")) %>%
                        dplyr::mutate(acc = ifelse(Individual_FramesNum == Correct_FramesNum, 1,0)) %>% 
                        dplyr::mutate(across(c("Trial", "Frame_Quantity"),.fns = as.numeric),
                                      acc = as.numeric(as.character(acc))) %>% 
                        dplyr::mutate(ResponseId = factor(ResponseId)) %>% 
                        dplyr::select(ResponseId,Trial,Frame_Quantity,Order,acc)%>% 
                        dplyr::mutate(task = "ID1")) %>% 
  dplyr::mutate(Order = as.numeric(Order)) %>% 
  na.omit()

mod1 <- glmer(acc ~ Order + (1|Trial)+(1|ResponseId), data = dfaccID1, family = binomial(link = "logit"))
summary(mod1)

dfsum <- summarySEwithin(data = dfaccID1, measurevar = "acc",withinvars = c("Order"), idvar= "ResponseId") %>% 
  arrange(Order) %>% 
  mutate(Order = as.numeric(as.character(Order)))

maxplot <- ggplot(dfsum, aes(x=Order, y=acc))+ #setting up x and y axes
  stat_eye(aes(y = acc,  dist = "norm", arg1 = acc, arg2 = se), show_interval = FALSE)+
  ylab("Accuracy")+ #naming y-axis
  xlab("Order of Frames")+ #naming x-axis
  ggtitle("Experiment 1 (Identify Max)")+  #naming the figure (or the plot)
  theme_tidybayes() + #theme of the plot
  scale_x_continuous(breaks=seq(2, 70, 4))+
  scale_y_continuous(limits=c(0,1.1),breaks=seq(0, 1, .1)) +
  geom_smooth(method = lm,se=F)+
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size=7, color =  "black"))
maxplot

dfsum <- summarySE(data = dfaccID1, measurevar = "acc",
                   groupvars = c("Order","ResponseId")) %>% 
  arrange(Order) %>% 
  mutate(Order = as.numeric(as.character(Order)))

maxplot <- ggplot(dfsum, aes(x=Order, y=acc))+ #setting up x and y axes
  stat_eye(aes(y = acc,  dist = "norm", arg1 = acc, arg2 = se), show_interval = FALSE)+
  ylab("Accuracy")+ #naming y-axis
  xlab("Order of Frames")+ #naming x-axis
  ggtitle("Experiment 1 (Identify Max)")+  #naming the figure (or the plot)
  theme_tidybayes() + #theme of the plot
  scale_x_continuous(breaks=seq(2, 70, 4))+
  scale_y_continuous(limits=c(0,1.1),breaks=seq(0, 1, .1)) +
  geom_smooth(method = lm,se=F)+
  facet_wrap(~ResponseId)+
  theme(plot.title = element_text(size = 12, hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 1, size=7, color =  "black"))
maxplot

