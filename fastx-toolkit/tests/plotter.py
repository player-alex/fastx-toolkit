# POWERED BY CHATGPT

import matplotlib.pyplot as plt
import numpy as np

# 데이터
labels = ['1M', '2.5M', '5M', '10M', '30M', '50M', '100M']
old_values = [5259, 12974, 27106, 52969, 157420, 262089, 528681]
buffering_values = [865, 2123, 4528, 8256, 24489, 40974, 85367]
openmp_values = [533, 1247, 2370, 4601, 13411, 21974, 52954]

# 바 그래프 위치 설정
x = np.arange(len(labels))

# 바 크기 설정
width = 0.3

# 그래프 그리기
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

# 각 데이터 추가 (시인성이 좋은 색상 선택)
bars1 = ax.bar(x - width, old_values, width, label='Old', color='#1f77b4')  # 파란색
bars2 = ax.bar(x, buffering_values, width, label='Buffering', color='#ff7f0e')  # 주황색
bars3 = ax.bar(x + width, openmp_values, width, label='OpenMP', color='#2ca02c')  # 초록색

# 바 위에 수치 추가
def add_value_labels(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 포인트 위에 텍스트 배치
                    textcoords="offset points",
                    ha='center', va='bottom')

# 각 바에 수치 추가
add_value_labels(bars1)
add_value_labels(bars2)
add_value_labels(bars3)

# 레이블 및 제목 설정
ax.set_ylabel('Time (ms)')
ax.set_title('FASTX Statistics: Performance Comparison')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

# 레이아웃 조정
plt.tight_layout()

# 그래프 표시
plt.show()
