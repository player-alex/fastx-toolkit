# POWERED BY CHATGPT

import matplotlib.pyplot as plt
import numpy as np

# ������
labels = ['1M', '2.5M', '5M', '10M', '30M', '50M', '100M']
old_values = [5259, 12974, 27106, 52969, 157420, 262089, 528681]
buffering_values = [865, 2123, 4528, 8256, 24489, 40974, 85367]
openmp_values = [533, 1247, 2370, 4601, 13411, 21974, 52954]

# �� �׷��� ��ġ ����
x = np.arange(len(labels))

# �� ũ�� ����
width = 0.3

# �׷��� �׸���
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

# �� ������ �߰� (���μ��� ���� ���� ����)
bars1 = ax.bar(x - width, old_values, width, label='Old', color='#1f77b4')  # �Ķ���
bars2 = ax.bar(x, buffering_values, width, label='Buffering', color='#ff7f0e')  # ��Ȳ��
bars3 = ax.bar(x + width, openmp_values, width, label='OpenMP', color='#2ca02c')  # �ʷϻ�

# �� ���� ��ġ �߰�
def add_value_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 ����Ʈ ���� �ؽ�Ʈ ��ġ
                    textcoords="offset points",
                    ha='center', va='bottom')

# �� �ٿ� ��ġ �߰�
add_value_labels(bars1)
add_value_labels(bars2)
add_value_labels(bars3)

# ���̺� �� ���� ����
ax.set_ylabel('Time (ms)')
ax.set_title('FASTX Statistics: Performance Comparison')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

# ���̾ƿ� ����
plt.tight_layout()

# �׷��� ǥ��
plt.show()
