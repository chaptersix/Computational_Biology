import datetime
# Dynamic Programming
#  - general and power algorithmic design technique
#  - Good for optimization problems
#  - careful bruteforce
#  - Richard Bellman in 1953 (Bellman Equation)


def NaiveFib(x):
    # correct but not efficient
    # lots of repeated work
    # Exponential Complexity!

    def fib(n):
        if n <= 2:
            f = 1
        else:
            f = fib(n-1) + fib(n-2)
        return f
    return(fib(x))


# Memoization (not memorization): to remember solutions already computed
def memoFib(x):
    # linear complexity
    # no repeat work
    # high space complexity
    # top down approach
    memo = {}

    def fib(n):
        if n in memo:
            return memo[n]
        if n <= 2:
            f = 1
        else:
            f = fib(n-1) + fib(n-2)
        memo[n] = f
        return f
    return(fib(x))

# Dynamic Programming is bottom up


def bottomUpFib(x):

    def fib(n):
        m = {}
        for k in range(1, n+1):
            if k <= 2:  # base case
                f = 1
            else:
                f = m[k-1] + m[k-2]
            m[k] = f
        return(m[n])
    return(fib(x))


num = 800

t = datetime.datetime.now()
a = memoFib(num)
t = datetime.datetime.now() - t
print("Memo execution time fib({}): {}\nAnswer: {}".format(num, t.microseconds, a))

t = datetime.datetime.now()
a = bottomUpFib(num)
t = datetime.datetime.now() - t
print("Bottoem-Up execution time fib({}): {}\nAnswer: {}".format(num, t.microseconds, a))

num = 40

t = datetime.datetime.now()
a = NaiveFib(num)
t = datetime.datetime.now() - t
print("Naive execution time fib({}): {}\nAnswer: {}".format(num, t.microseconds, a))
