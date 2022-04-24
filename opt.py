import subprocess
import optuna

def marathon(trial):
    c0 = trial.suggest_uniform('c0', 1, 1000)
    cmd = ["./marathon_main.exe", str(c0)]
    y = -float(subprocess.run(cmd, stdout=subprocess.PIPE).stdout) # 最小化しかないので、最大化したいときはマイナスをつける。
    print(y)
    return y

if __name__ == '__main__':
    study = optuna.Study(study_name='optuna', storage='sqlite:///optuna_database.db') # 最適化をいったんやめたり、並列化したいならこちら。
#   study = optuna.create_study() # 最適化を単発でやるだけなら、こちら
    study.optimize(marathon, n_trials=500)
