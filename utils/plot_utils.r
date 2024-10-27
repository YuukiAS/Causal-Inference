make_histogram = function(D, prop, type = "stacked") {
    #' Make a histogram of propensity score, showing the distribution of D=0 and D=1

    df = data.frame(as.factor(D), prop)
    if (type == "stacked" {
        ggplot(df, aes(x = prop)) +
            geom_histogram(data = subset(df, D == 0), aes(y = ..count.., fill = "D=0"), bins = 30, alpha = 0.8, color = "brown", position = "identity") +
            geom_histogram(data = subset(df, D == 1), aes(y = -..count.., fill = "D=1"), bins = 30, alpha = 0.8, color = "blue", position = "identity") +
            scale_y_continuous(labels = abs) + 
            labs(x = "Propensity Score", y = "Frequency") +
            theme_minimal() +
            scale_fill_manual(labels = c("D=0", "D=1"), values = c("brown", "blue")) + 
            guides(fill = guide_legend(title = NULL))  # remove title of legend
    }) else {
        # todo
        ggplot(data.frame(D=as.factor(D), prop), aes(x = prop, fill = D)) +
            # side-by-side histogram
            geom_histogram(position = "dodge", bins = 20, alpha = 0.8 , color = "white", size=3) +
            scale_fill_manual(values = c("blue", "red"),  labels=c("D=0", "D=1")) +
            theme_minimal() +
            labs(x = "Propensity Score", y = "Frequency", fill = "D")
    }
    
}