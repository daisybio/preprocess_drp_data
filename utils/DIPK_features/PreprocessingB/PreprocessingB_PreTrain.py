# code adapted from https://github.com/reinej03/DIPK_Bachelor_Thesis/tree/main and https://github.com/user15632/DIPK
import torch.optim as optim
from torch.utils.data import DataLoader
import joblib
import pandas as pd

from PreprocessingB_Model_DAE import *

# set parameters
seed = 14946934
lr = 1e-4
batch_size = 1018
EPOCHS = 1000
noising = True
save_path_log = 'PreTrain.txt'
save_path_model = 'PreTrain.pkl'
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

RMA_df = pd.read_csv('../Data/RMA.csv', sep='\t')
RMA_df_no_index = RMA_df.iloc[:, 1:]
GEF = torch.tensor(RMA_df_no_index.values, dtype=torch.float32)
# setup seed
sampler = setup_seed(seed)
# create model
encoder = Encoder(len(GEF[0])).to(DEVICE)
decoder = Decoder(len(GEF[0])).to(DEVICE)
loss_func = nn.MSELoss()
params = [
    {'params': encoder.parameters()},
    {'params': decoder.parameters()}
]
optimizer = optim.Adam(params, lr=lr)
# load data
my_collate = CollateFn()
train_loader = DataLoader(MyDataSet(GEF), batch_size=batch_size, shuffle=True, collate_fn=my_collate)
test_loader = DataLoader(MyDataSet(GEF), batch_size=batch_size, shuffle=False, collate_fn=my_collate)
# train model
for epoch in range(EPOCHS):
    # training
    encoder.train()
    decoder.train()
    epoch_loss = 0
    for it, Ft in enumerate(train_loader):
        Ft = Ft.to(DEVICE)
        if noising:
            z = Ft.clone()
            y = np.random.binomial(1, 0.2, (z.shape[0], z.shape[1]))
            z[np.array(y, dtype=bool)] = 0
            Ft.requires_grad_(True)
            output = decoder(encoder(z))
        else:
            output = decoder(encoder(Ft))
        loss = loss_func(output, Ft)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        epoch_loss += loss.detach().item()
    epoch_loss /= (it + 1)
    if epoch % 10 == 9:
        print('Epoch {}, loss {:.8f}'.format(epoch, epoch_loss))
        with open(save_path_log, 'a') as file0:
            print('Epoch {}, loss {:.8f}'.format(epoch, epoch_loss), file=file0)
    # saving
    if epoch % 1000 == 999:
        joblib.dump((encoder, decoder), save_path_model)
