# NetPro

Network Embedding framework in Multiplex Networks.

We consider node feature vector should contains all of information(Chemical,ATC,Target...)

The application of NetPro:

Drug combination prediction. Data from DrugBank, DrugCombDB, SIDER et al.

Five steps for NetPro:

1.Constructing the multiplex drug-drug association network

2.Performancing Random Walk with Restart to identify closely related drug pairs

3.Establishing randomforest models to predict synergistic drug combinations

4.Evaluating the toxicity risk for drug combinations with strong synergism

5.Predicting the therapeutic potential of drug combinations for different cancer types.
