# NetPro

Network propagation-based machine learning framework.

We consider node feature vector should contains all of information(Chemical,ATC,Target...)

The application of NetPro:

Drug combination prediction.  Data from DrugBank, DrugCombDB, SIDER et al.

Five steps for NetPro:

1. Constructing the multiplex drug-drug association network.

2. Performancing random walk with restart to identify closely related drug pairs.

3. Establishing random  forest models to predict synergistic drug combinations.

4. Evaluating the toxicity risk of drug combinations with strong synergy.

5. Predicting the therapeutic potential of drug combinations for different cancer types.
