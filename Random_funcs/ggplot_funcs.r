


# Rotate text
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#fix scale
p+scale_y_continuous(breaks = round(seq(min(test$n_eGenes_fdr_0_0_5), max(test$n_eGenes_fdr_0_0_5),by=300),1))


#position of labels

theme(axis.title.x = element_text(margin=margin(t=10)),axis.title.y = element_text(margin=margin(r=10)))


#basic barplot
p<-ggplot(data=df, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity")
p




### multi line


ggplot(data=res2, aes(x=PCs)) +
geom_line(aes(y=Astro,color="Astro"))+
geom_line(aes(y=Endo,color="Endo"))+
geom_line(aes(y=Excitatory,color="Excitatory"))+
geom_line(aes(y=Inhibitory,color="Inhibitory"))+
geom_line(aes(y=Microglia,color="Microglia"))+
geom_line(aes(y=Oligo,color="Oligo"))+
geom_line(aes(y=OPC,color="OPC"))+
geom_line(aes(y=Per,color="Per"))+
geom_line(aes(y=total,color="total"))+
geom_line(aes(y=unique,color="unique"))+
geom_point(aes(y=Astro),color=color_plane_1[1])+
geom_point(aes(y=Endo),color=color_plane_1[2])+
geom_point(aes(y=Excitatory),color=color_plane_1[3])+
geom_point(aes(y=Inhibitory),color=color_plane_1[4])+
geom_point(aes(y=Microglia),color=color_plane_1[5])+
geom_point(aes(y=Oligo),color=color_plane_1[6])+
geom_point(aes(y=OPC),color=color_plane_1[7])+
geom_point(aes(y=Per),color=color_plane_1[8])+
geom_point(aes(y=total),color=color_plane_1[9])+
geom_point(aes(y=unique),color=color_plane_1[10])+
scale_colour_manual(values=c("Astro"=color_plane_1[1],
"Endo"=color_plane_1[2],
"Excitatory"=color_plane_1[3],
"Inhibitory"=color_plane_1[4],
"Microglia"=color_plane_1[5],
"Oligo"=color_plane_1[6],
"OPC"=color_plane_1[7],
"Per"=color_plane_1[8],
"total"=color_plane_1[9],
"unique"=color_plane_1[10]))
