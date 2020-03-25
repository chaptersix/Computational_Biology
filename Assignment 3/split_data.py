import glob
import random
files = glob.glob('fasta/*')
# first 112 is the training data
files = [x.replace("fasta\\", "").replace(".fasta", "") for x in files]
random.shuffle(files)

training_proteins = files[:112]
testing_proteins = files[112:]

print(training_proteins)
#['2cua', '1nps', '2arc', '1jvw', '1im5', '1jbe', '1ctf', '1a3a', '1h0p', '1h2e', '1gbs', '1xkr', '1i58', '1whi', '1avs', '1svy', '1jwq', '1rw7', '1d0q', '1nb9', '1ny1', '1dlw', '1xdz', '1w0h', '2tps', '1jo0', '1bkr', '1chd', '1ek0', '1vjk', '1cjw', '1g2r', '1tqh', '1pko', '1mug', '1dqg', '1tqg', '1pch', '1hxn', '1g9o', '1lpy', '1jfx', '1rw1', '1i1j', '1c52', '1fcy', '1d1q', '1m4j', '1nrv', '1xff', '1dbx', '1wkc', '1j3a', '1fk5', '1tif', '1vmb', '1i4j', '1fl0', '1i1n', '1beb', '1jbk', '1jfu', '1h98', '1iwd', '1lm4', '1vfy', '1p90', '1a6m', '1gzc', '1hh8', '1iib', '1ej8', '2hs1', '1czn', '1fvk', '1k7c', '1qjp', '1ku3', '2vxn', '1aoe', '1jyh', '1fqt', '1jos', '5ptp', '1i5g', '1atl', '1r26', '1ne2', '1k7j', '1a70', '1hfc', '1f6b', '1cc8', '1atz', '1d4o', '1hdo', '1gmi', '1lo7', '1k6k', '1gmx', '1kid', '1fna', '1i71', '1dsx', '1guu', '1ql0', '1roa', '1ryb', '1cxy', '1t8k', '2mhr', '1eaz']
print("|||||||||")

print(testing_proteins)
# ['1vp6', '1fvg', '1dmg', '1jkx', '1bdo', '1ag6', '1gz2', '3bor', '1m8a', '1brf', '1ej0', '1aba', '1o1z', '1kw4', '1ktg', '3dqg', '1ihz', '1cke', '1jl1', '1dix', '1wjx', '2phy', '1jo8', '1kqr', '1mk0', '1vhu', '1fx2', '1htw', '1beh', '1aap', '1kq6', '1h4x', '1c9o', '1c44', '1tzv', '1qf9', '1bsg', '1smx'
