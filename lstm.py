import keras
import numpy as np
import os
import random
import tensorflow as tf
from collections import Counter

from keras.models import Sequential
from keras.optimizers import Adam
from keras.layers import Dense
from keras.layers import LSTM

model = Sequential()
dropout_frac = 0.0
learn_rate = 0.01

model.add(LSTM(units=100, dropout=dropout_frac, activation='tanh', recurrent_dropout=dropout_frac, return_sequences=True, input_shape=(600, 1)))
model.add(LSTM(units=100, dropout=dropout_frac, activation='tanh', recurrent_dropout=dropout_frac))
model.add(Dense(1, activation='sigmoid'))

adamopt = Adam(learning_rate=learn_rate)
model.compile(loss='binary_crossentropy', optimizer=adamopt, metrics=['accuracy'])

print('\nTrain...\n')
print(model.summary())

outdir = 'data'
nepoch = 4

val_x = np.fromfile(f'all.validation.50k.x.fixed.dat').reshape(-1,600,1)
with open(f'all.validation.50k.y.fixed.txt') as f:
    val_y = np.array([int(y.replace('.0','')) for y in f.read().splitlines()]).reshape(-1,1)

checkpoint_dir = 'rep_model_checkpoints'

for manual_epoch in list(range(nepoch)):
    print(f'epoch: {manual_epoch}')
    epoch_checkpoint_fp = f'./{checkpoint_dir}/dropout{dropout_frac}.epoch{manual_epoch}.last.keras'
    if os.path.exists(epoch_checkpoint_fp):
        model = tf.keras.models.load_model(epoch_checkpoint_fp)
        print(f'loading epoch_checkpoint_fp: {epoch_checkpoint_fp}')
        continue
        
    for chunk in range(20):
        print(f'chunk: {chunk}')
        checkpoint_fp = f'./{checkpoint_dir}/dropout{dropout_frac}.epoch{manual_epoch}.chunk{chunk}.keras'
        if os.path.exists(checkpoint_fp):
            model = tf.keras.models.load_model(checkpoint_fp)
            print(f'loading checkpoint_fp: {checkpoint_fp}')
            continue

        train_x = np.fromfile(f'{outdir}/chunk{chunk}.x.dat').reshape(-1,600,1)

        with open(f'{outdir}/chunk{chunk}.y.txt') as f:
            train_y = np.array([int(y) for y in f.read().splitlines()]).reshape(-1,1)
        
        #train_y = all_y
        if chunk == 0:
            print(np.shape(train_x))
            print(np.shape(train_y))
        #callback_list = [keras.callbacks.ModelCheckpoint(
        #    filepath=f'./checkpoints/model.{manual_epoch:02d}-{val_loss:.4f}', 
        #    save_freq='epoch', verbose=1, monitor='val_loss', 
        #    save_weights_only=True, save_best_only=False)]
        #model.fit(train_x, train_y, batch_size=500, epochs=1, callbacks=callback_list)
        model.fit(train_x, train_y, batch_size=200, epochs=1, validation_split=.2)
        model.save(checkpoint_fp)
        print(f'lastest model saved to: {checkpoint_fp}') #models were failing before the first epoch
    model.save(epoch_checkpoint_fp)
    print(f'end of epoch model saved to: {epoch_checkpoint_fp}')